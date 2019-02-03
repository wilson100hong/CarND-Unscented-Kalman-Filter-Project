# Unscented Kalman Filter Project Writeup
Self-Driving Car Engineer Nanodegree Program

---

[//]: # (Image References)

[image1]: ./images/data1_all.png "Dataset 1 - All Sensors"
[image2]: ./images/data1_laser.png "Dataset 1 - Laser"
[image3]: ./images/data1_radar.png "Dataset 1 - Radar"
[image4]: ./images/data2.png "Dataset 2 - All Sensors"
[image5]: ./images/nis_laser.png "Dataset 1 - NIS Laser"
[image6]: ./images/nis_radar.png "Dataset 1 - NIS Radar"

---

## Starter code
* [Udacity repo] (https://github.com/udacity/CarND-Unscented-Kalman-Filter-Project)

## Other Important Dependencies
* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF`

## Code Style

Please stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html) as much as possible.

## Project Rubrics

### Accuracy
1. px, py, vx, vy output coordinates must have an RMSE <= [.09, .10, .40, .30] when using the file: "obj_pose-laser-radar-synthetic-input.txt" which is the same data file the simulator uses for Dataset 1.

* The RMSE of dataset #1 at step 499 is in range

    | Variable | RMSE   |
    |----------|--------|
    | px       | 0.0692 |
    | py       | 0.0805 |
    | vx       | 0.3248 |
    | vy       | 0.2383 |

![Screenshot][image1]


* The RMSE of dataset #2 at step 498 has `vx` out of range though

    | Variable | RMSE   |
    |----------|--------|
    | px       | 0.0690 |
    | py       | 0.0683 |
    | vx       | 0.5625 |
    | vy       | 0.2089 |

![Screenshot][image4]


### Follows the Correct Algorithm
1. Your Sensor Fusion algorithm follows the general processing flow as taught in the preceding lessons.
Yes it follows the steps introduced from the course: predict then update measurement.

2. Your Kalman Filter algorithm handles the first measurements appropriately.
Yes in `ukf.cpp`, if the filter is not initialized, it will only update the state `x_` and `time_us`,
according to the types of sensor.

3. Your Kalman Filter algorithm first predicts then updates.
Yes in `ukf.cpp` function `ProcessMeasurement()`, it calls `Prediction()` and then `UpdateRadar()` or `UpdateLidar()`.

4. Your Kalman Filter can handle radar and lidar measurements.
Yes the filter will apply unscented transform to radar measurement. For lidar measurement, since the transform in linear,
I use same procedure as in normal KF. The experiment results shows same RMSE.


### Sensors Comparison
I compare the RMSE when only enabling radar or lidar. Both of them has worse RMSE than full UKF.

* The RMSE of dataset #1 when only use radar:

    | Variable | RMSE   |
    |----------|--------|
    | px       | 0.0692 |
    | py       | 0.0805 |
    | vx       | 0.3248 |
    | vy       | 0.2383 |

![Screenshot][image3]

* The RMSE of dataset #1 when only use lidar:

    | Variable | RMSE   |
    |----------|--------|
    | px       | 0.1506 |
    | py       | 0.2131 |
    | vx       | 0.3658 |
    | vy       | 0.2648 |

![Screenshot][image2]

### NIS Analytics
I tune the noises:

 | Noise     | value|
 |-----------|------|
 | std_a     | 1.5  |
 | std_yawdd | 0.5  |

The dataset #1 NIS values are under `data/`

1. Radar 

    | Percentile | Value  |
    |------------|--------|
    | 95%        | 7.2879 |
    | 5%         | 0.2682 |
 
![Screenshot][image6]
 
2. Laser
 
    | Percentile | Value  |
    |------------|--------|
    | 95%        | 4.6738 |
    | 5%         | 0.1319 |
 
![Screenshot][image5]


### Bonus challenge
See [https://github.com/wilson100hong/CarND-Catch-Run-Away-Car-UKF](https://github.com/wilson100hong/CarND-Catch-Run-Away-Car-UKF)