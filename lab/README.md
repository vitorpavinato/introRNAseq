# Introduction to RNA-seq

We are going to use R in our lab. For a very short introduction to R, take a look at this video from [DataCamp](https://www.youtube.com/watch?v=SWxoJqTqo08&ab_channel=DataCamp). If you are interested in dive into R a little more, I recommend going through the material in [this R crash course](https://billpetti.github.io/Crash_course_in_R/). Those materials will help you feel more comfortable with the code presented in the lab.

## Installation
For the lab, you need to install R and Rstudio. First install R (that is the language) and then install Rstudio (that is a environment, and IDE, to run R in a better interface). Follow the instructions for your computer provided in the links below.
- [R](https://cran.r-project.org); 
- [Rstudio](https://www.rstudio.com);

I am aware that you would need administrative access to install Rstudio in OSU's computers if you try to install manually (download from Rstudio site.)
However, you should be able to install if you do it through **Software Center**. Follow the instruction in [here](https://guides.osu.edu/tdai/intro) to install R and Rstudio OSU's computers. You can find the Software Center by typing it in the search bar at left bottom of Windows Desktop screen. I tried in one of our computer's lab and it worked. 

Also, install the following R packages before the lab (follow the instructions provided in the package page):
- [deseq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html);
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html);
- [gplots](https://www.rdocumentation.org/packages/gplots/versions/3.1.1).

An example how to install one of the packages (repeat the same steps for the other, making sure you copy the code from each of above links) 
1. Open Rstudio (make sure you have R and Rstudio properly installed in your computer).
2. Go to the command prompt (the one with ">" o the screen);
3. Paste the command found in the link. 

For example, for DESeq2 package, run these commands:
```dotnetcli
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
```

And to install edgeR package, run these commands:
```dotnetcli
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
``` 

4. Hit enter. You should see some prompts on your screen;
5. Repeat the same for `edgeR` and `limma`; 
6. For `gplots`, find the prompt and type:
```dotnetcli
install.packages("gplots")
``` 
7. Hit enter.

If everything went well, you are ready for the lab. 