---
title: "Random Proposal Code"
output: html_notebook
---


Create random data(Well distributed) to show why we need our variables to
normalized
```{r}
ex_data  = rnorm(n = 100, mean = 0, sd = 1)
```

```{r}
hist(ex_data)
```

Compare with population data (Maybe even orignal data)

```{r}
df_og = read.csv("C:\\Users/Laura/OneDrive - Western Sydney University/Spring 2022/Mathematics Project/Data Sets Modified/2015_done.csv")

df_pop = read.csv("C:\\Users/Laura/OneDrive - Western Sydney University/Spring 2022/Mathematics Project/Data Sets Modified/Population_done.csv")
```

```{r}
head(df_og, 5)
```


```{r}
head(df_pop,5)
```


```{r}
hist(df_og$dist_from_syd)
```
Highly skewered.
