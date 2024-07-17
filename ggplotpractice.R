# load libraries
install.packages("ggplot2")
library("ggplot2")

data(iris)
IrisPlot <- ggplot(iris, aes(Sepal.Length, Petal.Length, colour=Species)) + geom_point()
print(IrisPlot + labs(x = "Sepal length (cm)", y = "Petal length (cm)") + ggtitle("Petal and sepal length of iris"))

# add text
IrisPlot + annotate("text", x = 6, y = 5, label = "text")

# add repeat
IrisPlot + annotate("text", x = 4:6, y = 5:7, label = "text")

# highlight an area
IrisPlot + annotate("rect", xmin = 5, xmax = 7, ymin = 4, ymax = 6, alpha = .5)

# segment
IrisPlot + annotate("segment", x = 5, xend = 7, y = 4, yend = 5, colour = "black")

#---------------------------------------

bp <- ggplot(PlantGrowth, aes(x=weight, y=group)) + geom_point()

ggplot(iris, aes(Sepal.Length, Petal.Length, colour=Species)) + geom_point(shape=1) + geom_smooth(method=lm, se=FALSE)


# MACHINE LEARNING WORKFLOW

# EXPLORING AND PREPROCESSING----------
# Check for missing values/variables with little to no variation
# Exploratory analysis: Graphs & Summary Statistics
# outlier analysis
# normalization and standardization
# upsampling and downsampling easiest ways to combat an unbalanced dataset

# ASSESSING ML MODELS------------------
# Metric: MSE or RMSE -> Mean Square Error
# Metrics for assessing classification model: Accuracy, Sensitivity or True Positive Rate, Specificity or True Negative Rate
# Overfitting: problem, we should estimate how well the model performs out of sample
# Training and test data set: training contains 50% or more of all the available data
# cross-validation: available data set is randomly divided into k groups of ~equal size. all k groups used as test set. can be repeated t times
# LOOCV: Leave One Out Cross Validation: test set contains one data point, model trained on all data points except the one left out
# BOOTSTRAPPING: estimate a population parameter, estimates out-of-sample error. bootstrapped samples are used as training data

# LINEAR REGRESSION MODELS-------------

# Decision trees are the most basic models for classification problems
