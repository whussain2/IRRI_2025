###############################################################################
# Course: Fundamentals of Genomic Prediction-2025
# Module 1: Introduction and Learning R. 
# Authors: Waseem Hussain and  Mahender Anumalla 

###############################################################################


# Installing Packages

  install.packages("ggplot2")
  install.packages("AGHmatrix")

#Install these Packages
# AGHmatrix
# BGLR
# rrBLUP
# sommer
# lme4
# rmarkdown

## # Installs BiocManager if not installed
##   if (!require("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
##   BiocManager::install(version = "3.20")
## # Now Install R Package
##   BiocManager::install("GWASTools")

## # First Install Devtools
##   install.packages("devtools")
## # Then Intsall Package tidyr
##   devtools::install_github("tidyr")

# Examples of R Expression
# Add
  2+2
# Subtract
  4-2
# Square
    2^2
# Division
    10/2
# Multiple
    5*5
# Log
    log(10)
# Exponential
    exp(10)

# Create sample data
  x <- seq(1, 10, by=1)  # x values from 1 to 10
  y <- x^2  # y values are the squares of x

# Plot the data
  plot(x, y, 
       main="Plot of y = x^2",  # Title of the plot
       xlab="X-axis",            # Label for the x-axis
       ylab="Y-axis",            # Label for the y-axis
       col="blue",               # Color of the points
       pch=19)   
# Question: How many Arguments are in function plot()?
# Can you figure positional matching

# RATIONAL OPERATORS (Returns TRUE or FALSE)
# Is 3 > 4
   3>4
# Is 5 > 4
  5>4
# Is 5==5
  5==5
# Is 2 not equal to 3
  2!=3
# LOGICAL OPERATORS
  vector <- c(1, 2, 3, 4, 5) # Creating a simple vector
# Check which elements are greater than 2 and less than 5
  result <- vector[vector > 2 & vector < 5]  # result will be c(3, 4)
  result

# NUMERIC TYPE
  num1 <- 5     # Integer
  num2 <- 3.14   # Decimal (floating-point)
  num2
# INTEGER TYPE
  int_num <- 10L # Integer type
  int_num
# CHARACTER TYPE
  char<- "Hello, R!"  # Character string
  char
# LOGICAL
  True_type <- TRUE
  False_type <- FALSE
  False_type
# COMPLEX 
complex_num <- 3 + 2i  # 3 is the real part, and 2 is the imaginary part
complex_num
#Assignment operator in R is <-, why not Equal sign =
#Question: Could you use class() function to check the type of data

# Examples of numeric scalers
    a <- 100
    a
    b <- 3 / 100
    b
    c <- (a + b) / b
    c
# Examples of character scalers
    d <- "ship"
    d
    e <- "cannon"
    e
    f <- "Do any modern armies still use cannons?"
    f
# Could you add two scalers of type character?

# Example of Vector
  X<-c(1,2,3,4,5, 6) # five components
  X
# Check length
  length(X)
# Creating a vector of 2s eight times
    d<-rep(2, 8)
    d
    d<-rep(c(1,2,3,4), 5)
    d
    d<-rep(c(1:4), 2)
    d
#INDEXING or SUBSETTING
    
    x[1] # Extract first element
    x[c(1,2, 6)] # Extract first, second and 6 element in vector x (sub vector)
    x[1:3] # Extracts elemnts from 1 to thrid position 
    x[-6] # Drops last element
    sum(x) # sum function adds all elements
    mean(x) # get mean of x
    min(x) # get minimal value
    max(x) # Get maximum value

# Creat a matrix
## 2.3.1 First approach
    m2<-matrix(1:9, nrow = 3, ncol = 3) # Matrix with 3 rows and columns
    m2 # Matrix filled by rows
    m2<-matrix(1:9, nrow=3, ncol=3, byrow = TRUE)
# Assign names to rows and columns
    dimnames(m2)<-list(c("X","Y","Z"), c("A","B","C"))
    m2
# Access and chnage the row and column names
    colnames(m2) # Get column names
    row.names(m2)
# Change the namaes of columns
    colnames(m2)<-c("A.1", "B.1", "C.1")
    m2
# Change the name of just first column
    colnames(m2)[1]<-"B.B"
    m2
## Second approach to create the matrix

# We can use cbind() and rbind() functions, column and row bind
    x<-c(4,2,3,6) # create a vector of x
    y<-c(8,5,6,9) # create a vector of y
    m.col<-cbind(x,y) # now use cbind to bind two vectors column wise
    m.col
    m.row<-rbind(x,y) # now use cbind to bind two vectors row wise
    m.row
### Extract the elements of matrices. 
## Square bracket [] indexing method. Elements can be accessed as var[row, column].
# First create a new matrix
    m2<-matrix(1:12, nrow=3, ncol=4, byrow = TRUE)
    class(m2)
    m2[1,] # Extracts first row 
    m2[,4] # Extracts  4th column
    m2[1,4] # Extracts first elemnt in row 1 and column 1
    m2[c(1,3), c(1,2)]
    m2[-1,] # Leaves first row
## Modify the matrix
    m2[1,1]<-10 # Chnages single element in first row and first column
    m2
    m2[m2>11]<-20 # Chnage the elements in matrix greater than 12
    m2
# Add colum or row to existing matrix
    m2<-cbind(m2, c(4,8,12))
    m2
    x<-c(4,2,3,6,7) # create a vector of x
    m2<-rbind(m2,x)
    m2
# Note dimensions should be same to add vectors
    m2<-rbind(m2, c(1,2))
    m2
# Check the dimensions of matrix
    dim(m2)
# Transpose the matrix
    t(m2)

# Create a data.frame
    Genotypes <- c("Genotyp1", "Genotype1", "Genotype2", "Genotype2")
    Replication <- c("1","1","2","2")
    Block<- c("Block1", "Block2", "Block1", "Block2")
    Yield<-c(2500,3500,3200,4500)
    mydata <- data.frame(Genotypes,Replication,Block,Yield)
    class(mydata)
# Check the structure
    str(mydata)
# Change the varaibales
    mydata$Block<-as.factor(mydata$Block)
    mydata$Replication<-as.factor(mydata$Replication)
# Check the structure again
    str(mydata)
    levels(mydata$Replication) # Determine the number of levels for replication.
# Check column names
    names(mydata)
# Subsetting or Indexing
    mydata[1:2,]
    mydata[c(1,3),]
# Select columns using $ sign
    mydata$Yield # select Yield column
    mydata[mydata$Yield>3500,] # Extract row that has yield greater than 3500
# Subset based on factor levels
    help("subset")
    mydata2<-subset(mydata,Block=="Block1") # select just Block1
    mydata2
    mydata2<-subset(mydata,Block=="Block1" & Replication=="1") # select just Block1
    mydata2

# Creating a simple list
  my_list <- list(name = "Waseem", age = 39, height = 5.5, is_student = FALSE)
  my_list
# Mixed List: Creating a list with different types of elements
  mixed_list <- list(vector = c(1, 2, 3), 
                   matrix = matrix(1:4, nrow = 2), 
                   char = "Hello", 
                   logical = c(TRUE, FALSE))
  mixed_list
# Creating a nested list
  nested_list <- list(person1 = list(name = "Alice", age = 25), 
                    person2 = list(name = "Bob", age = 30))
# Extracting using $
# Accessing elements in the nested list
  person1_name <- nested_list$person1$name
  person2_age <- nested_list[["person2"]][["age"]]

# Example 1
    x <- 10
    y <- 12
    if(x < y) {
      print("x is less than 12!")
    }
# Example 2
    data("iris")
    if(mean(iris$Sepal.Length)!=mean(iris$Petal.Length)){
      print("It is true")
    }

# Example 1
    x <- 10
    y <- 12
    if(x > y) {
      print("x is greater than y")
    } else {
      print("y is greater than x")
    }
# Example 2
    if(mean(iris$Sepal.Length)<mean(iris$Petal.Length)) {
      print("This is true and mean of sepal length is less than petal length")
    } else {
      print("This is false and mean of petal length is less than sepal length")
    }

# Example 1
    x <- 12
    y <- 13
    if(x > y) {
      print("x is greater")
    } else if(x < y) {
      print("y is greater")
    } else {
      print("x and y are equal")
    }
# Example 2
    x <- 12
    y <- 12
    if(x > y) {
      print("x is greater")
    } else if(x < y) {
      print("y is greater")
    } else {
      print("x and y are equal")
    }

# Example of Multiple conditions
# First let us check means
    mean(iris$Sepal.Length)
    mean(iris$Petal.Length)
    mean(iris$Sepal.Width)
    mean(iris$Petal.Width)
# Test the multiple conditions using operators
    if(mean(iris$Sepal.Length)>mean(iris$Petal.Length) &&
       mean(iris$Sepal.Width)>mean(iris$Petal.Width) ) {
      plot(iris$Sepal.Length, main="Sepal length",
           ylab="Values")
    } else {
      plot(iris$Petal.Length, main="Petal length",
           ylab="Values")
    }

# IFELSE Examlple
# Example for ifelse function
    iris$mycolumn<-ifelse(iris$Sepal.Length>iris$Petal.Length, 
                          "TRUE", "FALSE") 
    iris$mycolumn<-ifelse(iris$Sepal.Length<iris$Petal.Length,
                          iris$Sepal.Length/iris$Petal.Length, NA)  
# Let us change the values of sepal length
    iris$Sepal.Length<-iris$Sepal.Length/3
    iris$mycolumn2<-ifelse(iris$Sepal.Length>iris$Petal.Length,
                           iris$Sepal.Length/iris$Petal.Length, NA)  
### Apply funtion
    
# First create a matrix
    my.matrx <- matrix(c(1:10, 11:20, 21:30), nrow = 10, ncol = 3)
    my.matrx
    colnames(my.matrx)<-c("Length", "Breadth", "Width")
#Get sum across rows by using apply function
    sumrow<-apply(my.matrx, 2, mean)
    sumrow
# Creating own function
    sumrow<-apply(my.matrx, 1, function (x) sum(x)*2)
    sumrow
# Creating function outside and then apply
# Creating a function to calculate Cofficient of variation
    my.cofvar<- function(x){
      (sd(x)/mean(x))*100
    }
# Now apply it it dataframe or matrix
    cv<-apply(my.matrx,1, my.cofvar)
    cv

## lapply Function
# Create a list first
    list.1<-list(Length=c(6,4,8,6.5),breadth=c(7,8,6,8),width=c(6.8,7.2,6.6,8))
    list.1
# Get mean of all lists using lapply function
    mean.all <- lapply(list.1,mean)
    mean.all # Mean for all the lists

## sapply Function
# Create a vector first
    vec<-c(1,2,3,4,5,6,7,8)
    mean.vec <- sapply(vec,mean)
    mean.vec
# Additional argument simplify
    vec<-c(1,2,3,4,5,6,7,8)
    mean.vec <- sapply(vec,mean, simplify = FALSE) # RETURNS REULTS AS LAPPLY
    mean.vec

# tapply Function
# First creating a simple data frame
    sample.data <- data.frame(Genotypes=c("Geno1","Geno2","Geno3","Geno4",
                                          "Geno5","Geno6","Geno7"),
    Yield=c(240, 220, 211, 230, 203, 241, 212),
    Type=factor(c("LR","LR","LR","CV", "CV","CV","LR")))
# Now apply tapply function to get mean for various factor levels
    mean<-tapply(sample.data$Yield, sample.data$Type, sum)
    mean

#  For loop
# Example 1
    for (i in 1:5) { 
      print(i^3) 
    }
# Example 2
    x <- c(-8, 9, 11, 45) 
      for (i in x) { 
        y<-x/2
        print(y) 
      } 
# Example 3
# Histogram for iris data all columns
# Get iris data
    data(iris)
# Histograms
    par(mfrow = c(2, 2)) # Create 2 x 2 plotting matrix
    for (i in 1:ncol(iris[,c(1:4)])){ # loop over columns
      plot <- hist(iris[,i], main=paste("Trait", i))
      plot
    }  
# Your Assignment to run a simpl example of WHILE AND REPEAT LOOPS

# Load library
   library(stringr)
# First create a string
    x<-c("ICAR", "IRRI", "Training")
# Get number of characters
    nchar(x)
# Convert to lower case
    y<-tolower(x)
    y
# Convert to upper case
    y<-toupper(x)
    y
# Upper or lower case conversion with casefold()
# upper case folding
    y<-casefold(x, upper = TRUE)
    y
# Replace 'IRRI' by 'IRRIGO'
    help("chartr")
    z<-"IRRI ICAR TRAINING"
    z2<-chartr("I", "G",z)
    z2
# abbreviate species name column in iris data set
    iris$Species<- abbreviate(iris$Species, minlength = 3)
    head (iris$Species)
#with package stringer
    head(abbreviate(iris$Species, minlength = 4, dot = FALSE, 
                     strict = FALSE))
# Replace substrings with substr()
    help("substr")
    substr(iris$Species, 1, 3)<-"MY"
# Check similar function in stringr
# Concatenating strings    
x<-c("ICAR", "Colloboration", "Genomic")
y<-c("IRR", "Training", "Prediction")
yx<- paste(x, y)
yx
# Without space
yx<- paste0(x, y)  # Output: "HelloWorld"
yx

# First create a two character vectors
    set1 = c("ICAR", "IRRI", "TRAINING", "PROGRAM")
    set2 = c("IRRI", "TRAINING", "PROGRAM", "2020")
# union of set1 and set2
    union(set1, set2) #discards any duplicated values in the provided vectors
#Set intersection with intersect()
    intersect(set1, set2) # common between two
# Set difference with setdiff()
    setdiff(set2, set1)
# First create set 3 and set 4
    set3=c("IRRI", "ICAR", "TRAINING", "PROGRAM")
# Set equality with setequal()
    setequal(set1, set3)
# Exact equality with identical()
    identical(set1, set3)
# Element contained with is.element()
# Create elements
    elem1 = "IRRI"
    elem2 = "2020"
# elem1 in set 1
    is.element(elem1, set1) #if an element is contained in a given set of character strings 
# elem1 in set10?
    elem1 %in% set1 # alternative
# Sorting with sort()
#sort() allows us to sort the elements of a vector, either in increasing order (by
#default) or in decreasing order using the argument decreasing:
# sort (decreasing order)
    set5<-c("boy", "girl", "apple")
    set6<-c("boy", "girl", "apple", "2020")
    sort(set5)
    sort(set6)
#Repetition with rep()
# repeat 'x' 4 times
    help("rep")
    rep("Apple", 4)
    help("paste")
# Check similar  functions in stringr

# Replacing first occurrence with sub()
# string
  help("sub")
  Rstring = c("The r Foundation",
            "for Statistical Computing",
            "R is FREE software",
            "r is a collaborative project")
# substitute 'R' with 'R
  sub("R", "RR", Rstring, ignore.case = TRUE)
#Replacing all occurrences with gsub()
# substitute
  gsub("R", "RR", Rstring)

## Splitting Character Vectors
  sentence = c("R is a collaborative project with many contributors")
# split into words
  strsplit(sentence, " ")
# telephone numbers
  tels<-c("5_1_0_548_20")
  tels = c("510-548-2238", "707-231-2440", "650-752-1300")
# split each number into its portions
  strsplit(tels, "_")
# Check similar  functions in stringr: your assignment

## Additional on Regular expressions
# Meta-characters
# string
  money = "$money"
# the naive but wrong way
  sub(pattern = "$", replacement = "", x = money)
# the right way in R
  sub(pattern = "\\$", replacement = "_", x = money)
# Sequences
#replaces the first match, while gsub() replaces all the matches.
# replace digit with '_'
  sub("\\d", "_", "the dandelion war 2010")

  gsub("\\d", "_", "the dandelion war 2010")

  sub("\\D", "_", "the dandelion war 2010")

  gsub("\\D", "_", "the dandelion war 2010")

#Spaces and non-spaces
# replace space with '_'
  sub("\\s", "_", "the dandelion war 2010")

  gsub("\\s", "_", "the dandelion war 2010")

# replace non-space with '_'
  sub("\\S", "_", "the dandelion war 2010")

  gsub("\\S", "_", "the dandelion war 2010")

# Setting a new working directory
  setwd("~/Documents/Research/Workshops/IRRI-IIRR_2024/Module1")
  #setwd("C:Documents/Research/Workshops/") # In windows, replace path
# Get a working directory
  getwd()
# Create a new working directory
  dir.create("Test")

# Import .Txt file
# Read .csv file data into R. File in in Data folder
  mydata<-read.csv("./Data/iris.txt", header = TRUE)
#Check this function to know more help("read.table")
# Export the file as .csv file
  write.csv(mydata, "./Data/iris_export.csv",
            row.names=TRUE)

# Import .Txt file
# Read tabular data into R
  mydata<-read.table("./Data/iris.txt", header = TRUE, sep = "\t")
#Check this function to know more help("read.table")
# Export the file as .txt file
  write.table(mydata, "./Data/iris_export.txt",
              append=FALSE,sep="/t",row.names=TRUE,col.names=TRUE)

# IMPORT EXCEL SHEETS
# Load the Library
  library(readxl)
# Read the excel file
# Specify sheet by INDEX
  my_data <- read_excel("./Data/iris.xlsx", sheet = 1)
# Specify sheet by Name
  my_data <- read_excel("./Data/iris.xlsx", sheet = "IRIS_COPY")
# EXPORT EXCEL FILES
  library(writexl)
# Writing the Existing file
  #write_xlsx(my_data, "./Data/myIrisExoprt.xlsx")
# Check with Library ("xlsx")

# Read from web directly
  rice.pheno <- read.table("http://www.ricediversity.org/data/sets/44kgwas/RiceDiversity_44K_Phenotypes_34traits_PLINK.txt", 
    header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# First five rows and columns
  rice.pheno[1:5, 1:5]

# Save an .rds object in folder
  saveRDS(my_data, file = "./Data/my_data.rds")
# Restore  or read the same object
  myrds<-readRDS(file = "./Data/my_data.rds")

# Install packages and load
  packages = c("dplyr", "tidyr", "reshape", "reshape2")
# Create a function which will install the package if it is not installed
  package.check <- lapply(packages, FUN = function(x) { # apply lapply to list of packages
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE) # install dependencies if required'library(x, character.only = TRUE) # Load the package
      }
    })

# Read iris data file
  iris<-read.csv(file="./Data/iris.csv", header = TRUE)
  iris$Species.abr<-abbreviate(iris$Species, minlength = 3)
######################PACKAGEDPLYR#############################################  

# Function FILTER   
  iris1<-iris %>% filter(Species == "setosa", Species.abr == "sts")
# same in base R
  iris2<-iris[iris$Species == "setosa" & iris$Species.abr == "sts", ]
# Function ARRANGE   
  iris1<- iris1 %>% arrange(Sepal.Length)
# Use desc() to order a column in descending order:
    iris1<-iris1%>% arrange(desc(Sepal.Length))
#Choose rows using their position with slice() 
    iris %>% slice(1:5)
# in base R
    test<-iris[c(1,5, 6:10), ]
# Drop
    iris1<-iris %>% slice(-c(1:5))
# Function SELECT
# Select columns by name
    iris1<-iris %>% select(Sepal.Length, Sepal.Width)
# Select all columns between 
    iris1<-iris %>% select(Sepal.Length:Petal.Width)
# Select all columns except one
    iris1<- iris %>% select(!Sepal.Length)
    #iris %>% select(!c(Sepal.Length, Petal.Width))
# Select all columns ending with length
    iris1<-iris %>% select(starts_with("Species"))
    iris1<-iris %>% select(ends_with("abr"))
# more: starts_with(), ends_with(), matches() and contains() (Your assignment)
    #select(iris, contains("Length"))
    iris1<-iris %>% select(contains("."))
    #iris1<-select(iris, contains("."))
# Function MUTATE   
# Add single column    
    iris1<-iris %>% mutate(Additional = Sepal.Length/ 100)
# Add several columns
    iris1<-iris%>%mutate(Additional2 = sqrt(Sepal.Length), Check = "vrs")
# What is alternative in R base function (Your Assignment)
# Function SUMMARISE 
  summary<-iris %>% summarise(Sepal.Length = mean(Sepal.Length, na.rm = TRUE))
  summary
  iris %>% summarise(Mean=mean(Sepal.Length, na.rm = TRUE), SD=sd(Sepal.Length, na.rm = TRUE))
# Summaries with column names
  iris %>% summarise(mean(Sepal.Length, na.rm = TRUE),mean(Petal.Length, na.rm = TRUE))
# Summarise a subset of rows
  iris%>%
    slice(1:10) %>%
    summarise(sum(Sepal.Length))
#Function DISTINCT
  iris%>% distinct(across(Species))
  
###############################PIPE MULTIPLE TASKS#############################
# Pipe to perform multiple tasks
  test<-iris %>%
    filter(Species == "setosa")%>% # filter
    mutate(SePAL.2=Sepal.Length/ 10)%>% # create new column
    select(SePAL.2, Species) # and select the column
# Another pipe example
    test2<- iris%>%                   
    group_by(Species) %>%    # Group by these variables
    summarise(  # Summarise data
    length.mean = mean(Sepal.Length),      
    width.mean = median(Sepal.Width),
    n = n()                    
    ) 
# Assignmet: Check with Tidyr, reshape2, data.table functions

#SCATTER PLOT
  plot(iris$Sepal.Length,
# Add features
  xlab="Genotype", # Add x label name
  ylab="Sepal length",  # add y lable name
  main="Scatter plot for Sepal length",# add main title
  type="b",
# Different types of graphsx
  las=1, #0 is the default, with text always parallel to its axis.
  pch=20, 
  cex=1, 
  col="blue",
  lty=2, # line type
  lwd=4, # line or point width
  font.lab = 1, # font type
  col.lab="black", 
  cex.lab = 1.5, # size of axis
  xlim=c(4,8), ylim=c(1,7), # adjust limits of axis
)

# BAR PLOTS
  par(mfrow = c(1,2)) # Plot in One row with two Columns
  barplot(iris$Sepal.Length, 
  main="Bar Plot for Sepal length",
  xlab="Sepal length", ylab="Frequency",
  col="blue", 
  width = 5,
  horiz = FALSE)
  #?barplot
# do more on bar plots
  barplot(cbind(Sepal.Length, Petal.Length)~Lines, data=iris, 
        subset = Species=="setosa",
        col = c("lightblue", "mistyrose"),
        main="Bar Plot for Sepal length",
        xlab="Sepal length", ylab="Frequency")

  par(mfrow = c(1,2)) 
  hist(iris$Petal.Length, 
     freq=TRUE,
     breaks=12,
     col="red",
     xlab="Sepal length",
     main="Histogram for Seapl length")
  rug(jitter(iris$Petal.Length)) #   add rugs
  lines(density(iris$Petal.Length), col="blue", lwd=2) # add density curve
# Another graph
   hist(iris$Sepal.Length, 
     freq=FALSE,
     breaks=8,
     col="blue",
     xlab="Sepal length",
     main="Histogram for Seapl length")
  rug(jitter(iris$Sepal.Length)) #   add rugs
  lines(density(iris$Sepal.Length), col="red", lwd=4) # add density curve

#BOX PLOTS  
boxplot(Sepal.Length ~ Species, data=iris,
        main="Box plot for Sepal length",col="blue",
        xlab="Species",
        ylab="Values")


# SCATTER PLOT
# Install and Load the libraries
  packages = c("ggplot2", "grid","ggthemes","plotly")
# Create a function which will install the package if it is not installed
  package.check <- lapply(packages, FUN = function(x) { # apply lapply to list of package 
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE) # install dependencies if required
      library(x, character.only = TRUE) # Load the package
      }
    })
# Use ggplot Function to plot
  ggplot(mtcars, aes(x=wt, y=mpg))+ 
    geom_point()+
    geom_point(size=2, shape=23, color="red")+ # Change the point size, and shape
    #geom_smooth(method=lm)+  # Add the regression line
    #geom_smooth(method=lm, se=FALSE)+ # Remove the confidence interval
    geom_smooth(method=lm, se=TRUE, linetype="dashed",
              color="darkred") # Change the line type and color
# SCATTER PLOT MULTIPLE GROUPS
# Change point shapes by the levels of cyl
  mtcars$cyl<-as.factor(mtcars$cyl)
  ggplot(mtcars, aes(x=wt, y=mpg, shape=cyl, color=cyl, size=cyl)) +
    geom_point()+
    geom_smooth(method=lm)+
    geom_smooth(method=lm, se=FALSE, fullrange=FALSE)
# Change point shapes and colors manually
  ggplot(mtcars, aes(x=wt, y=mpg, color=cyl, shape=cyl))+
    geom_point()+
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
    scale_shape_manual(values=c(3, 16, 17))+ 
    scale_color_manual(values=c('blue','red', 'green'))+
    theme(legend.position="top")+
    geom_rug()+
    theme_few() 
# Change the point sizes manually
  ggplot(mtcars, aes(x=wt, y=mpg, color=cyl, shape=cyl))+
    geom_point(aes(size=cyl)) + 
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
    scale_shape_manual(values=c(3, 16, 17))+ 
    scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))+
    scale_size_manual(values=c(2,3,4))+
    theme(legend.position="top")

  #png(file="myBar_plot.png", width=14, height =10, 
      #units = 'in',res=700)  # function to Export as .png
  #tiff(""myBar_plot.tiff", units="in", width=5, height=5, res=500) # export as tiff
 # pdf(""myBar_plot.pdf", width = 4, height = 4) # Export as PDF
# ggplot 2 function
  ggplot(mydata, aes(x = Species, y = Sepal.Length, fill=Species))+
    #geom_bar(stat = "identity", width=0.5, color="black")+
    geom_bar(stat = "identity", width=0.5, 
             color="black", position=position_dodge())+
    theme_few()+ # Select the theme
    ylim(0, 10)+ # Set Y limit
    labs(title = "", x = "Species", y = "Sepal.length")+ # add titles
    theme (axis.title.x = element_text(color="blue", size=18), # modify axis titles
           axis.title.y = element_text(color="blue", size=18))+
    theme(axis.text= element_text(face = "bold", color = "darkred", size =18))+
    scale_fill_manual(values=c("#CC79A7", "#009E73", "#e79f00"))+ # Manual coloring
    theme(legend.title = element_text(colour="black", size=12),
          legend.position = "right", 
          legend.text = element_text(colour="red", size=14))
  #dev.off()

# Hiostograms in ggplot2
ggplot(mydata, aes(Sepal.Length)) +
  geom_histogram(color="darkblue", fill="lightblue", bins = 10)+
  theme_light()+
# adjust x values and breaks
  #geom_histogram(breaks=seq(4, 8, by =0.5), color="darkblue", fill="lightblue")+
  geom_vline(aes(xintercept=mean(Sepal.Length)), # adding the line to represent mean
             color="darkred", linetype="dashed", size=1)+
  labs(title="",x="Value", y = "Count")+
  theme (plot.title = element_text(color="black", size=14, face="bold", hjust=0),
         axis.title.x = element_text(color="black", size=10, face="bold"),
         axis.title.y = element_text(color="black", size=10, face="bold")) +
  theme(axis.text= element_text(face = "bold", color = "black", size = 8))


# Create Your Function
# Build function
my.cofvar<- function(x){
  (sd(x)/mean(x))*100
}
# Create a data frame
sample.data2 <- data.frame(Genotypes=c("Geno1","Geno2","Geno3","Geno4",
                                       "Geno5","Geno6","Geno7"),
                           Yield=c(240, 220, 211, 230, 203, 241, 212),
                           Height=c(110, 140, 160, 135, 125, 145, 150),
                           Type=factor(c("LR","LR","LR","CV", 
                                         "CV","CV","LR")))

# Now calculate CV for Yield and Height variable
my.cofvar(sample.data2$Yield)
my.cofvar(sample.data2$Height)


################################END############################################