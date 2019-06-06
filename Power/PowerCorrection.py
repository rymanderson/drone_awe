class PowerCorrection:
        'Class used to make corrections to the drone\'s power requirement based on maneuvers, weather, payload, etc.'

        # class variables go here:
        addedpower      = 0.0           # watts

        # methods go here:
        def __init__(self, sink):
            # sink is an object describing the reason for an increase/decrease in the power requirement (e.g., an instance of the Rain.py class)
                

print("Successfully imported `PowerCorrection.py`")