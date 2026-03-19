/* rw_cascade_simple.c
   Cascaded reaction-wheel control using IMU angle & rate directly.
   - Outer PID:  angle θ -> rate setpoint ω_ref   (slow)
   - Inner PI:  rate ω  -> motor command u       (fast)
*/

#include <stdint.h>
#include <math.h>
#include <LSM6DS3.h>
#include <Wire.h>
#include <SPI.h>
#include <nRF24L01.h>
#include <RF24.h>

//Create an instance of class LSM6DS3 (IMU)
LSM6DS3 myIMU(I2C_MODE, 0x6A);    //I2C device address 0x6A
//Create an instance of class RF24 (Radio nRF24L01)
RF24 radio(3, 4); // CE, CSN 

/* ====== Config ====== */
#define M_PI 3.14159265358979323846
#define PREC             5.0          // Precision +/- degrees
#define T_THRES          1.0          // Thruster kick-in angle rate threshold [deg/s]
#define T_IMP            100          // Time of thruster impulse [ms] 
#define T_CORR           1.5          // Time correction factor of IMU when thrusters kick in
#define CTRL_HZ          1000.0f      // inner loop frequency
#define OUTER_HZ          100.0f      // outer loop frequency (CTRL_HZ divisible by OUTER_HZ)
#define TS               (1.0f/CTRL_HZ)
#define OUTER_DIV        ((int)(CTRL_HZ/OUTER_HZ))

/* Limits */
#define PWM_MAX            0.8f*255       // |u| <= PWM_MAX (Assuming saturation occurs beyond 80% duty cycle)
#define RATE_REF_MAX       10             // [rad/s]
#define RATE_REF_SLEW      20             // [rad/s^2] limit on ω_ref change

/* Inner (rate) PI gains — tune on your rig */
#define KR_P               58.8f 
#define KR_I               4.0f 

/* Outer (position) runtime-tunable PID gains — much slower than inner */
volatile float Kp_pos   = 0.8f;
volatile float Ki_pos   = 0.05f;
volatile float Kd_pos   = 0.10f;
volatile float Kd_alpha = 0.10f;  // for filtered D term (0..1)

const byte address[6] = "00001";

#pragma pack(push,1)
struct CmdAck {
  float theta_deg;
  float kp;
  float ki;
  float kd;
  float kd_alpha;
  uint8_t seq;
};
#pragma pack(pop)

/* ====== Helpers ====== */
static inline float clampf(float x, float lo, float hi){ return x<lo?lo:(x>hi?hi:x); }
static inline float wrap_pi(float a){	// round off angle θ to within 2*PI rad (360 degrees) 
  while(a >  M_PI) a -= 2.0f*M_PI;
  while(a < -M_PI) a += 2.0f*M_PI;
  return a;
}

/* ====== Controller state ====== */
static float theta_ref = 0.0f;   // commanded angle [rad]
static float omega_ref = 0.0f;   // commanded body rate [rad/s]
static float omega_ref_prev = 0.0f;

static float integ_pos = 0.0f;   // outer loop integrator
static float integ_rate = 0.0f;  // inner loop integrator
static float omega_d_filt = 0.0f;
static int   outer_div = OUTER_DIV;

int flag = 0;
long timeStart, timeEnd, dTime;
float aX, aY, aZ, gX, gY, gZ;
float gZ1;
float dgZ = 0.0;
float ang = 0.0;
float preAng = 0.0;
float pregZ = 0.0;
float bias = 0; 
float offset = 0;
float biasgZ = 0; 
const float accelerationThreshold = 0.0; // threshold of significant in G's
const int numSamples = 119;
int samplesRead = numSamples;

const uint8_t N = 3;                 // how many float-point variables (measurements) we send
float payload[N];

struct IMUMeasurements {
  float ang;  // angle (deg)
  float gZ;   // gyro Z (deg/s)
};


void setup() {
  // put your setup code here, to run once:
  integ_pos = integ_rate = 0.0f;
  outer_div = OUTER_DIV;
  theta_ref = 0.0*M_PI/180.0; //-90: left; 90: right
  omega_ref = omega_ref_prev = 0.0f;

  // put your setup code here, to run once:
  pinMode(0, OUTPUT);                  // IN1 pin for RW motor
  pinMode(1, OUTPUT);                  // IN2 pin for RW motor
  pinMode(D5, OUTPUT);                  // ENA pin for thruster 1 (Left Thruster)
  pinMode(D6, OUTPUT);                  // ENA pin for thruster 2 (Right Thruster)

  digitalWrite(0, LOW);                // Turn off RW motor initially
  digitalWrite(1, LOW);
  digitalWrite(D5, LOW);                // Turn off thrusters initially
  digitalWrite(D6, LOW);

  myIMU.begin();

  if (!myIMU.begin())
  {
    biasgZ = -myIMU.readFloatGyroZ()/5000;	// Bias for gyro reset at start 
  }

  ang = 0;
  preAng = ang;
  gZ= 0;
  pregZ = gZ;

  radio.begin();
  radio.openWritingPipe(address);
  radio.setPALevel(RF24_PA_MIN);

  radio.setAutoAck(true);
  radio.enableAckPayload();
  radio.setRetries(5, 15);

  radio.stopListening();
}

void loop() {
  // put your main code here, to run repeatedly:
  /* Measurements from IMU (already in angle & rate) */
  //timeStart = millis();
   const uint32_t t0 = micros();
  // call the IMU function
  IMUMeasurements meas = readIMU();

  float theta = wrap_pi(meas.ang * M_PI / 180.0f);
  float omega = meas.gZ   * M_PI / 180.0f;
  float e_theta;

  /* ---- Outer PID: θ -> ω_ref (runs at OUTER_HZ) ---- */
  if(--outer_div <= 0){
    outer_div = OUTER_DIV;

    e_theta = wrap_pi(theta_ref - theta);
    float up = Kp_pos * e_theta;

    // simple conditional anti-windup: only integrate when not saturating at the rate limit
    float ui = integ_pos + Ki_pos * (1.0f/OUTER_HZ) * e_theta;
    // D term on measurement to avoid derivative kick; filter gyro rate
    omega_d_filt += Kd_alpha * (omega - omega_d_filt);
    float ud = -Kd_pos * omega_d_filt;
    float v  = up + ui + ud;

    // limit rate reference and apply slew-rate limit
    float v_sat = clampf(v, -RATE_REF_MAX, RATE_REF_MAX);
    float max_step = RATE_REF_SLEW / OUTER_HZ;
    float step = clampf(v_sat - omega_ref_prev, -max_step, max_step);
    omega_ref = omega_ref_prev + step;

    // update integrator only if not clamped hard
    if (v == v_sat) integ_pos = ui;  // accept integration when not saturating

    omega_ref_prev = omega_ref;
  }

  /* ---- Inner PI: ω -> u (runs at CTRL_HZ) ---- */
  float e_omega = omega_ref - omega;
  float up_rate = KR_P * e_omega;

  // conditional integration to avoid windup against PWM limits
  float ui_rate = integ_rate + KR_I * TS * e_omega;
  float u = up_rate + ui_rate;
  float u_sat = clampf(u, -PWM_MAX, PWM_MAX);
  if (u_sat > 0)
  {
    digitalWrite(0, LOW);
    digitalWrite(1, HIGH);
  }
  else if (u_sat < 0)
  {
    digitalWrite(0, HIGH);
    digitalWrite(1, LOW);
  }
  else
  {
    digitalWrite(0, LOW);
    digitalWrite(1, LOW);
  }
  if (u == u_sat) integ_rate = ui_rate;  // integrate only when not saturating

  /* Command for rocket (fan) motors during reaction wheel saturation */
  analogWrite(2, max(abs((int)u_sat), 0));
   
  const uint32_t t1 = micros();
  dTime = (t1 - t0) / 500.0f;   // ms for your integrator
  preAng = meas.ang;
  pregZ  = meas.gZ;

  if (flag == 0 && u_sat <= -PWM_MAX+1 && e_theta < -PREC*M_PI/180.0 && pregZ < -T_THRES) 
  {
    digitalWrite(D5, LOW);      //HIGH
    digitalWrite(D6, HIGH);
    delay(T_IMP);
    digitalWrite(D5, LOW);      
    digitalWrite(D6, LOW);
    dTime = dTime + T_CORR*T_IMP;
    flag = 1;
  }
  if (flag == 0 && u_sat <= -PWM_MAX+1 && e_theta > PREC*M_PI/180.0 && pregZ > -T_THRES) 
  {
    digitalWrite(D5, HIGH);      //HIGH
    digitalWrite(D6, LOW);
    delay(T_IMP);
    digitalWrite(D5, LOW);      
    digitalWrite(D6, LOW);
    dTime = dTime + T_CORR*T_IMP;
    flag = 1;
  }
  if (flag == 0 && u_sat >= PWM_MAX-1 && e_theta > PREC*M_PI/180.0 && pregZ > T_THRES)  
  {
    digitalWrite(D5, HIGH);      //LOW
    digitalWrite(D6, LOW);
    delay(T_IMP);
    digitalWrite(D5, LOW);      
    digitalWrite(D6, LOW);
    dTime = dTime + T_CORR*T_IMP;
    flag = 1;
  }
  if (flag == 0 && u_sat >= PWM_MAX-1 && e_theta < -PREC*M_PI/180.0 && pregZ < T_THRES)  
  {
    digitalWrite(D5, LOW);      //LOW
    digitalWrite(D6, HIGH);
    delay(T_IMP);
    digitalWrite(D5, LOW);      
    digitalWrite(D6, LOW);
    dTime = dTime + T_CORR*T_IMP;
    flag = 1;
  }
  if (fabs(e_theta) >= PREC*M_PI/180.0)
  {
    flag = 0;
  }

  /* Radio signal transmission */
  static uint32_t t_tx = 0;
  if (millis() - t_tx >= 50) {   // 20 Hz telemetry (adjust as needed)
    t_tx += 50;

    payload[0] = meas.ang;       // deg
    payload[1] = meas.gZ;        // deg/s
    payload[2] = (int)u_sat;   // PWM cmd (or normalize to [-1..1])

    bool ok = radio.write(payload, sizeof(payload));
    bool queued = radio.writeFast(payload, sizeof(payload)); // queue to FIFO
    radio.txStandBy(2);   // wait <= 2 ms max, then give up
    if (!queued) radio.flush_tx();

    // ---- If Disp queued an ACK payload, read θ_ref (in degrees) ----
  if (ok && radio.isAckPayloadAvailable()) {
    //float theta_deg_in;
    CmdAck incoming, last;
    // In case multiple ACK payloads are queued, drain them (keep the last)
    while (radio.isAckPayloadAvailable()) {
      radio.read(&incoming, sizeof(incoming));
      last = incoming;
    }
    // Convert degrees -> radians and wrap to [-pi, pi]
    theta_ref = wrap_pi(last.theta_deg * (M_PI / 180.0f));

    // 2) Update outer-loop gains with simple sanity clamps
    auto clampf = [](float x, float lo, float hi){ return x<lo?lo:(x>hi?hi:x); };
    Kp_pos   = clampf(last.kp,       0.0f, 10.0f);
    Ki_pos   = clampf(last.ki,       0.0f, 10.0f);
    Kd_pos   = clampf(last.kd,       0.0f, 10.0f);
    Kd_alpha = clampf(last.kd_alpha, 0.0f,  1.0f);
  }
  }

}

struct IMUMeasurements readIMU() {
  IMUMeasurements m;
  
  gZ = myIMU.readFloatGyroZ() + biasgZ - 0.43; // Gyro (z) offset to eliminate gyro drifting
  m.gZ = gZ;

  // integrate to get angle
  m.ang = preAng + 0.5f * (pregZ + gZ*0.55 + bias) * dTime / 1000.0f;

  return m;
}