# Estimating prey activity curves using a quantitative model based on a priori distributions and predator detection data
Daniel J. Herrera, Daniel Levy, Austin M. Green, William F. Fagan


### Background
Ecologists routinely study the impact of predators on prey activity patterns through the largely qualitative approach of comparing overlapping activity density plots. While this method offers some insight into predator-prey dynamics, it does not actually quantify the degree to which a predator's activity impacts prey activity. Here, we present a novel model that overcomes this shortcoming by using predator detections and an ideal prey activity curve to quantify the impact of predator activity on prey activity patterns. The model assumes that species strive to adhere to an ideal activity distribution and quantifies the degree to which a predator prompts a departure from this ideal curve. This approach improves our ability to quantify and predict predator-prey interactions as they pertain to activity patterns. For a more detailed justification, please see our paper published in Ecological Modelling (2024).
<br>
### Model Structure
Our model is premised upon the idea that prey species try to adhere to a specific activity pattern in the absence of predators. Coincidentally, these patterns happen to resemble familiar data distributions. For instance, let's assume that our prey species' activity curve (blue) resembles the centered beta distribution (alpha = 2, beta = 2):
![readme plot1](https://github.com/Dan-Herrera/Estimating_Activity_Curves/assets/66024392/26cc4f9d-634c-42c7-aa45-5f7e03366f2c)

<br>
Now let's suppose a predator (red) is active during the same time period, but does not adhere to the same activity pattern as the prey.

![readme plot2](https://github.com/Dan-Herrera/Estimating_Activity_Curves/assets/66024392/c85f88a8-42b8-4e6e-af73-a1544921a4b8)

<br>
The previous graph assumes that the prey (blue) is not impacted by the predator (red). But our model uses the proportional intensity of predator activity to alter the prey activity curve (with a lower bound of zero since negative activity is not possible). Additionally, the model acknowledges that prey will respond differently to various predators. Thus, we estimate the sensitivity (s) of the prey to the predator. If the prey is highly sensitive, then its activity curve will be more dampened than if the prey has a low sensitivity value. Here's three different outcomes where s = 0.25 (top), 0.50 (middle), and 0.75 (bottom).

![readme plot3](https://github.com/Dan-Herrera/Estimating_Activity_Curves/assets/66024392/2d8dae5c-5ff0-4fce-9c56-3f7dccff6760)

<br>
Simply removing activity from the prey's curve is not biologically realistic, however, since the animal still needs to be out foraging, searching for mates, etc. Thus, the model allows the prey to respond by compensating for activity lost earlier in the active period. Importantly, the model assumes compensation for lost activity becomes more important later in the active period since there are progressively fewer opportunities to make up for lost activity as the day goes on (assuming being active outside of the active period is not possible). In the following plot, the dashed blue line represents the prey's activity when s = 0.25 and no lost activity is compensated for. The solid blue line represents the prey's activity when s = 0.25 and lost activity is compensated for. Note that compensation attempts at time t seek to restore the realized activity curve to the ideal activity curve as some later time t+x. In other words, allowing for activity compensation reveals an activity curve that seeks to return (or closely resemble) the ideal a priori activity curve.

![readme plot4](https://github.com/Dan-Herrera/Estimating_Activity_Curves/assets/66024392/8b843c5f-c1c2-4b70-8752-8bea5670c093)

<br>

### Running the model on real data

The above plots were made using fictitious data, which aren't particularly exciting or compelling. But our model allows for the estimation of s when using real data. To fit your own data to this model, be sure to download the functions script from the code folder and source the script so the functions exist in your environment.
```{}
source("./YOUR FILE PATH/activityFunctions.R")
```

<br>
You'll want to use your own data, but for the sake of example we'll use the data that we used in the paper. You can access the entire data set via [Cove, et al. 2019](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.3353), our you can access a miniature version of this data set using the `UtahData()` function.


```{}
UtahData() #loads data into environment and names it "UtahData"
```
<br>
Note that these observation times are already in *POSIXlt* format. Use the `strptime()` function to convert your own data's timestamp column to *POSIXlt* format.

<br>
To convert observation timestamps into activity curves, use the `format_observedPDF()` function. If your goal is to predict prey distribution based on a known s value, you only need to use this function on the predator data. If s is unknown and you wish to estimate it, use this function for both predator and prey data. Note that time.start and time.stop must be the same across all species and functions.

```{}
cottontail <- format_observedPDF(df = UtahData, #name of data frame
                                 sp.col = "species", #name of species column in data frame 
                                 sp.name = "Mountain Cottontail", #name of focal species
                                 time.start = "18:00", #start of active period as character without date
                                 time.stop = "12:00", #end of active period as character without date
                                 time.col = "obs.time") #name of column that holds observation time stamps

fox <- format_observedPDF(df = UtahData,
                                 sp.col = "species", 
                                 sp.name = "Red Fox",
                                 time.start = "18:00",
                                 time.stop = "12:00",
                                 time.col = "obs.time")
```
<br>
To set up the a priori distribution, use the `format_idealPDF()` function.

```{}
ideal.dist <- format_idealPDF(time.start = "18:00", #specify start of active period as character without date
                              time.stop = "12:00", #specify end of active period as character without date
                              observedPDF = cottontail, #specify observed activity distribution
                              dist = "triangular", #specify distribution
                              param1 = 0, #various shape parameters governing the distribution (error messages will help you through these)
                              param2 = 1,
                              param3 = 0.2)
```
<br>
For simplicity, we will only make one a priori distribution here. But be sure to test several options in your own analysis (and read the paper to understand why). The `format_idealPDF()` function offers a number of different distributions that will seamlessly work with the other functions included in this script. To see which pre-set distributions are possible, feed the function any value that you know isn't a distribution to throw a helpful error:

```{}
ideal.dist <- format_idealPDF(time.start = "18:00",
                              time.stop = "12:00",
                              observedPDF = cottontail,
                              dist = "i don't know")
```
> ERROR: the following shapes are permissable (but custom shapes can be used as long as the vector length matches that of the disturbance vector and the vector sums to 1): beta, binomial, bimodal, kumaraswamy, normal, triangular, uniform.
<br>
Now that we have our ideal a priori distribution, as well as our observed distributions for both predator and prey, we can estimate the sensitivity of the prey to the predator using the estimate_activity() function.

```{}
pred.list <- list() #create an empty list
pred.list[[1]] <- fox #fill the first slot of the list with fox data

results <- estimate_activity(idealPDF = ideal.dist, #specify ideal a priori distribution
                             disturbance = pred.list, #specify list of predator activity distributions (even if only one predator)
                             observedPDF = cottontail, #specify observed prey activity distribution
                             confInt = TRUE, #indicate that we want to calculate a confidence interval
                             nIter = 1000, #specify the number of iterations to bootstrap over (this can take a while, so be patient or pick a small number)
                             AIC = TRUE) #instruct R to calcualte AIC score (only useful if comparing models)
```
| S | CI95 | negLL |  AIC  |
| :---         |    :---        |         :---   |:---        |
| 0.09791847   | (0.0884 - 0.1145)     | -531.6364    |1081.273 |

<br>
While plotting activity curves can be useful, the above-printed results are the real pride of this model. Recall that s is the sensitivity of the prey to the predator. Rather than telling us the degree of overlap between the predator and prey activity curves, estimating s allows us to determine the degree to which the predator actually changes a prey's activity. In other words, this method provides us mechanistic results rather than descriptive results.

<br>

That being said, we do still want to visualize our results. You can feed the estimated s value into a predictive model using the `predict_activity()` function. Similar to how predator activity curves must be provided as a list, sensitivity values do as well.

![readme plot5](https://github.com/Dan-Herrera/Estimating_Activity_Curves/assets/66024392/c88b8d67-4a2c-4861-bc36-95f5436d2fd4)

<br>

This is slightly different than the output you'll see when using the functions on your own computer. Here, we see the ideal distribution (triangular distribution) in grey, the observed cottontail distribution in red, and the predicted cottontail distribution in blue. As you can see, the predicted line resembles the observed line, and sufficiently deviates from the a priori distribution. In other words, it works!

<br>

While this function will plot a rudamentary visualization in your plotting window, it will also save the results to your environment as a data frame (provided you run the function as a named object, eg. name <- function(arguments)). You can then use that data frame to create a more robust visualization, which is what I have shown here.

<br>

We hope that this model suits your needs, and that the functions included in the activityFunctions.R script are flexible enough to meet your needs. Feel free to use this model as is, feed the model different types of data (e.g., precipitation data instead of predator data), or improve upon the model by including additional terms. Also, if you're a better coder than I am, I would love help implementing a better likelihood estimator (the current estimator searches over values much faster than my homemade estimator could, but but crashes when I try to use two predator data sets in the same model)!

<br>

Finally, we encourage you to read the paper, and cite us if you decide on using this model. Please feel free to reach out if you have any questions or comments (herrerawildlife@gmail.com). You can also reach me via my [personal website](https://www.herrerawildlife.com/).
