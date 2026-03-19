"""
One-Axis Reaction Wheel Stick - MicroPython conversion

Converted from Arduino .ino to MicroPython for ESP32 / Raspberry Pi Pico.
Adjust pin assignments and PWM parameters to match your board and wiring.

Hardware:
  - MPU6050 IMU (I2C)
  - Nidec 24H brushed DC motor with driver
  - 3S LiPo battery with voltage divider on ADC pin
  - Active buzzer
"""

import struct
import sys
import math
import time
from machine import Pin, PWM, I2C, ADC, UART

# ---------- MPU6050 registers ----------
MPU6050_ADDR = 0x68
ACCEL_CONFIG = 0x1C
GYRO_CONFIG = 0x1B
PWR_MGMT_1 = 0x6B
PWR_MGMT_2 = 0x6C

# ---------- Pin assignments (adjust for your board) ----------
BRAKE_PIN = 8
PWM_PIN = 9
DIRECTION_PIN = 4
BUZZER_PIN = 12
VBAT_PIN = 35  # ADC-capable pin; adjust as needed

# ---------- PWM settings ----------
PWM_FREQUENCY = 20000  # 20 kHz

# ---------- Control gains ----------
X1 = 75.0   # Proportional (angle)
X2 = 5.25   # Derivative (angular velocity)
X3 = 0.04   # Integral (motor speed)
LOOP_TIME_MS = 10  # Control loop period in milliseconds

# ---------- Sensor config ----------
ACC_SENS = 0   # 0 = 2g, 1 = 4g, 2 = 8g, 3 = 16g
GYRO_SENS = 1  # 0 = 250 deg/s, 1 = 500 deg/s, 2 = 1000 deg/s, 3 = 2000 deg/s
GYRO_AMOUNT = 0.996  # Complementary filter gyro weight

# ---------- IMU offsets ----------
AcX_offset = -750
AcY_offset = 360
AcZ_offset = 0
GyZ_offset = 0

# ---------- Filter ----------
alpha = 0.40  # Low-pass filter coefficient for gyro

# ---------- State variables ----------
pwm_s = 0
motor_speed = 0
robot_angle = 0.0
Acc_angle = 0.0
gyroZfilt = 0.0
vertical = False

AcX = 0
AcY = 0
GyZ = 0

# ---------- Hardware init ----------
i2c = I2C(0, scl=Pin(22), sda=Pin(21), freq=400000)  # Adjust pins for your board

brake = Pin(BRAKE_PIN, Pin.OUT)
direction = Pin(DIRECTION_PIN, Pin.OUT)
buzzer = Pin(BUZZER_PIN, Pin.OUT)
vbat = ADC(Pin(VBAT_PIN))

motor_pwm = PWM(Pin(PWM_PIN))
motor_pwm.freq(PWM_FREQUENCY)

# PWM duty range: MicroPython uses 0-65535 (16-bit) or 0-1023 (10-bit)
# depending on port. We normalise to 0-65535.
PWM_MAX = 65535


def constrain(val, min_val, max_val):
    return max(min_val, min(max_val, val))


def map_value(x, in_min, in_max, out_min, out_max):
    return int((x - in_min) * (out_max - out_min) / (in_max - in_min) + out_min)


# ---------- I2C helpers ----------
def write_to(device, address, value):
    i2c.writeto_mem(device, address, bytes([value]))


def read_raw_word(device, reg):
    """Read a signed 16-bit big-endian value from two consecutive registers."""
    data = i2c.readfrom_mem(device, reg, 2)
    value = struct.unpack(">h", data)[0]
    return value


# ---------- IMU ----------
def angle_setup():
    global GyZ_offset

    write_to(MPU6050_ADDR, PWR_MGMT_1, 0)  # Wake up MPU6050
    time.sleep_ms(100)
    write_to(MPU6050_ADDR, ACCEL_CONFIG, ACC_SENS << 3)
    write_to(MPU6050_ADDR, GYRO_CONFIG, GYRO_SENS << 3)
    time.sleep_ms(100)

    # Calibrate gyro offset
    offset_sum = 0
    for _ in range(1024):
        angle_calc()
        offset_sum += GyZ
        time.sleep_ms(5)
    GyZ_offset = offset_sum >> 10

    # Beep twice to signal calibration complete
    buzzer.value(1)
    time.sleep_ms(70)
    buzzer.value(0)
    time.sleep_ms(80)
    buzzer.value(1)
    time.sleep_ms(70)
    buzzer.value(0)

    print("GyZ offset value =", GyZ_offset)


def angle_calc():
    global AcX, AcY, GyZ, robot_angle, Acc_angle, vertical

    # Read accelerometer X and Y
    AcX = read_raw_word(MPU6050_ADDR, 0x3B)
    AcY = read_raw_word(MPU6050_ADDR, 0x3D)

    # Read gyroscope Z
    GyZ = read_raw_word(MPU6050_ADDR, 0x47)

    # Apply offsets
    AcX += AcX_offset
    AcY += AcY_offset
    GyZ -= GyZ_offset

    # Integrate gyro to get angle
    robot_angle += GyZ * LOOP_TIME_MS / 1000.0 / 65.536

    # Accelerometer angle
    Acc_angle = math.atan2(AcY, -AcX) * 57.2958  # rad to deg

    # Complementary filter
    robot_angle = robot_angle * GYRO_AMOUNT + Acc_angle * (1.0 - GYRO_AMOUNT)

    # Vertical detection with hysteresis
    if abs(robot_angle) > 9:
        vertical = False
    if abs(robot_angle) < 0.3:
        vertical = True


# ---------- Motor ----------
def set_pwm(duty):
    """Set motor PWM duty cycle (0 to PWM_MAX)."""
    motor_pwm.duty_u16(constrain(int(duty), 0, PWM_MAX))


def motor_control(pwm_val):
    if pwm_val <= 0:
        direction.value(0)
        pwm_val = -pwm_val
    else:
        direction.value(1)
    # Map 0-255 to PWM_MAX-0 (inverted)
    set_pwm(map_value(pwm_val, 0, 255, PWM_MAX, 0))


# ---------- Battery ----------
def batt_voltage(voltage):
    if 8 < voltage <= 9.5:
        buzzer.value(1)
    else:
        buzzer.value(0)


# ---------- Serial tuning ----------
def tuning(uart):
    global X1, X2, X3

    if not uart.any():
        return
    time.sleep_ms(2)
    param = uart.read(1)
    if param is None:
        return
    param = chr(param[0])

    if not uart.any():
        return
    cmd = uart.read(1)
    if cmd is None:
        return
    cmd = chr(cmd[0])

    if param == 'p':
        X1 += 1 if cmd == '+' else (-1 if cmd == '-' else 0)
        print_values()
    elif param == 'i':
        X2 += 0.01 if cmd == '+' else (-0.01 if cmd == '-' else 0)
        print_values()
    elif param == 's':
        X3 += 0.005 if cmd == '+' else (-0.005 if cmd == '-' else 0)
        print_values()


def print_values():
    print("X1: {:.2f} X2: {:.2f} X3: {:.3f}".format(X1, X2, X3))


# ---------- Main ----------
def main():
    global pwm_s, motor_speed, gyroZfilt

    # UART for serial tuning (use UART0 or USB serial depending on board)
    uart = UART(0, baudrate=115200)

    # Initial PWM and pin state
    set_pwm(map_value(400, 0, 400, PWM_MAX, 0))
    buzzer.value(0)
    brake.value(1)
    time.sleep_ms(1000)

    print("PWM_FREQUENCY:", PWM_FREQUENCY)
    angle_setup()

    previous_t1 = time.ticks_ms()
    previous_t2 = time.ticks_ms()

    while True:
        current_t = time.ticks_ms()

        # Control loop (every LOOP_TIME_MS)
        if time.ticks_diff(current_t, previous_t1) >= LOOP_TIME_MS:
            tuning(uart)
            angle_calc()

            if vertical:
                brake.value(1)
                gyroZ = GyZ / 131.0  # Convert to deg/s
                gyroZfilt = alpha * gyroZ + (1 - alpha) * gyroZfilt
                pwm_s = -constrain(
                    int(X1 * robot_angle + X2 * gyroZfilt + X3 * -motor_speed),
                    -255, 255
                )
                motor_control(pwm_s)
                motor_speed += pwm_s
            else:
                motor_control(0)
                brake.value(0)
                motor_speed = 0

            previous_t1 = current_t

        # Battery check (every 500 ms)
        if time.ticks_diff(current_t, previous_t2) >= 500:
            raw_adc = vbat.read()  # 0-4095 on ESP32, 0-65535 on Pico
            batt_voltage(raw_adc / 53.7)
            previous_t2 = current_t


if __name__ == "__main__":
    main()
