# Development Lab Notebook

This document is meant to document design decisions made during the development of this package.

## Simple Power model 20190621

### Camera Power

In order to determine whether or not camera power should be accounted for, the power consumption of a GoPro was researched on https://cam-do.com/pages/power-consumption-by-gopro-camera-model. It was determined to be negligible compared to the rotor power, as the extreme case--a 4K camera in video capture--required only 5 watts of electricity. When added to the rotor power consumption, the effect on the curve was un-noticeable.

### Interpolate propulsive efficiencies

It was determined that the propulsive efficiency would change significantly with flight speed. The efficiency was calculated at two datapoints--one at the max endurance conditions and one at the max range conditions, as reported in the data sheet for the drone in question--based on momentum theory. These datapoints were then used to interpolate a new efficiency based on the velocity selected for the simulation according to: $$eta = eta1 - \frac{v1-v}{v1-v2} (eta1 - eta2)$$