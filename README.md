# MacArthur_model

## Matplotlib

Matplotlib is a Python library for creates plots.

## Pandas

Pandas is a fast, powerful, flexible and easy to use open source data analysis and manipulation tool, 
built on top of the Python programming language.

## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install matplotlib and pandas.

```bash
pip install matplotlib
pip install pandas 
```

### Usage

```python
import matplotlib.pyplot as plt
import pandas as pd

# Convert an array or list to Data Frame
saves_ni_m = pd.DataFrame(ni_m)

# Save the Data Frame as a csv file
saves_ni_m.to_csv('MacArthur Euler ni_mean '+'('+str(time)+').csv', index=False)

# Create a figure with two subplots
fig, axs = plt.subplots(1, 2)
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Authors and acknowledgment

Autor: Mauricio Silva Tovar

I would like to express my special thanks to my advisor Isaac Pérez Castillo who always encouraged me to give my best and supported me at all times.

## License

[MIT](https://choosealicense.com/licenses/mit/)
