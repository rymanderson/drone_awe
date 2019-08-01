# Ryan's Notebook

## How to predict lift and drag on the drone body

It has been observed in experimental data that power consumption decreases from the hover case at low velocities, and then increases at higher velocities. This phenomenon is to be modeled, but first must be somewhat understood.\\

First, it was hypothesized that lift on the drone airframe is causing this reduction in power. However, this notion is somewhat counter-intuitive since the airframe is at a negative angle of attack. Interestingly, this phenomenon was also observed in an IMAV paper on propellers without a drone body at varying thrusts and angle of attack. (http://www.imavs.org/papers/2017/321_imav2017_proceedings.pdf). They found that at low angles of attack, there is a reduction in power consumption at lower velocities compared to the hover case, while this phenomenon did not occur at high angles of attack. It is possible that the power reduction observed is more an artifact of propeller dynamics at varying angle of attack and velocity, and less to do with lift on the drone airframe.\\

To further investigate the mechanism associated with this effect, the drag and lift on a rotary drone is researched.

## Predicting rotor propulsion

Sources found include the following:

* https://vtol.org/files/dmfile/howHelosFly_presentation.pdf
* https://static1.squarespace.com/static/56e4a24bc6fc082c7577a416/t/5718f8ce45bf21d4d295ef24/1461254351404/Grant_AER1216_lecture.pdf

## Fixed wing aircraft

Sources may include:

* for maneuvers: https://pdfs.semanticscholar.org/579f/fb6a8688627f4da1a5380624fb630479bf25.pdf

## VTOL aircraft

Sources may include:

* design of a small eVTOL biplane: https://arxiv.org/pdf/1801.02938.pdf
