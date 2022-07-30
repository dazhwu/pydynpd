This document describes how to manipulate the outputs of the abond function. For instruction on what inputs are needed to run the abond function, please refer to [Tutorial](https://github.com/dazhwu/pydynpd/blob/main/vignettes/Tutorial.ipynb).

## Outputs of the abond function

After regression, the abond function returns a list of models:

![Figure 1](https://raw.githubusercontent.com/dazhwu/pydynpd/main/vignettes/Images/list_models.svg)

Multiple models could be returned because users may include "?" in command string. For example, in Example 8 in [Tutorial](https://github.com/dazhwu/pydynpd/blob/main/vignettes/Tutorial.ipynb), the following comand string generates 5 models:

    command_str='n L(1:?).n w k  | gmm(n, 2:3) pred(w k)'
    mydpd = regression.abond(command_str, df, ['id', 'year'])

### Usage:

regression.abond()

Quarto enables you to weave together content and executable code into a finished document. To learn more about Quarto see <https://quarto.org>.

## Running Code

When you click the **Render** button a document will be generated that includes both content and the output of embedded code. You can embed code like this:

```{python}
1 + 1
```

You can add options to executable code like this

```{python}
#| echo: false
2 * 2
```

The `echo: false` option disables the printing of code (only output is displayed).
