

## Models

Here are the independent models.


Each model is a class, and should be initialized (trained) with the test dataset upon initialization (`__init__` method).
Predictions are asked via its `Predict` method, taking a dataset row (FAT% results should be deleted from the row before prediction).
The preditction method should return a tuple of floats: (`FAT% Liver`, `FAT% Pancreas`).


By using different models we can evaluate how each approach performs, and take note
how modifications in models perform against its older versions and other models.
