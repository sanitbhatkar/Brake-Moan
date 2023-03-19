# Brake-Moan

## Motivation 
1. Brake moan noise is low-frequency noise occurring at low speed. Multiple factors like stick-slip, mode coupling and atmospheric conditions can be the cause. For the MotionSolve model, the input of the friction curve is required. 
2. Bajaj Auto uses three classes of Brake pads: a.Organic b. Sintered c. High COF Organic. The pads are procured from different vendors, thus the composition of the pad is different. Each pad has a different friction curve which affects the brake moan propensity due to stick slip. 
3. Experiments were conducted on the Brake Test Bench to understand the variation of the COF with temperature, pressure and disc velocity. The experiments were conducted for many brake cycles. Multi-class time series data is generated through the experiments.

## Current Work:
1. Pandas and NumPy Library is used for Data import, DataFrame storage, Aggregation and Data manipulation. Visualization is done with Matplotlib.
2. KNN is used to remove the outliers from the dataset. (PyOD library)
3. Regression of the data is done with the sci-kit learn library and some in-house codes. The model is trained on 60% data, tested on 20% data and validated with the remaining dataset for accuracy.
4. Dataset and the trained model is stored in pickle data format for easy access.

