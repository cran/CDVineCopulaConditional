

#' Random dataset from a given vine copula model
#' @description A random dataset simulated from a given 5-dimensional vine copula model. 
#' 
#'   
#' @format 
#' \describe{
#'   \item{\code{$data}}{An \code{1000 x 5} data set (format \code{data.frame}) with the uniform 
#'   variables \code{(U1,U2,U3,U4,U5)}.}
#'   \item{\code{$vine}}{\code{\link[VineCopula]{RVineMatrix}} object defyining the vine copula model from where \code{$data} was sampled.}
#' }
#' 
#' 
#' 
#' @examples
#' 
#' # Load data
#' data(dataset)
#' 
#' # Extract data
#' data <- dataset$data 
#' plot(data)
#' 
#' # Extract the RVineMatrix object from where the dataset was randomly sampled
#' vine <- dataset$vine
#' vine$Matrix 
#' vine$family
#' vine$par
#' vine$par2 
#' summary(vine)
#' 
#' 
#' 
#' @author Emanuele Bevacqua
"dataset"








#' Simulation from a conditional C- or D-vine
#' 
#' @description Simulates from a d-dimensional conditional C- or D-vine of the variables (\strong{Y},\strong{X}),  
#' given the fixed conditioning variables \strong{X}. The algorithm works 
#' for vines satysfying the requirements discussed in \emph{Bevacqua et al. (2017)}. The algorthm implemented here 
#' is a modified version of those form \emph{Aas et al. (2009)} and is shown in \emph{Bevacqua et al. (2017)}. 
#' 
#' 
#' @param RVM An \code{\link[VineCopula]{RVineMatrix}} object containing the information of the d-dimensional C- or D-Vine 
#' model (for further details about \code{\link[VineCopula]{RVineMatrix}} objects see the documentation file of the \code{VineCopula} package). 
#' If the full copula is 2-dimensional, RVM can be an \code{\link[VineCopula]{RVineMatrix}} object or a data.frame (or list) object 
#' where \code{$family}, \code{$par} and \code{$par2} are specified.
#' @param Condition A \code{N x Nx} matrix of the Nx conditioning variables.
#' For D-vine: data corresponding to the conditioning variable whose index is in \code{RVM$Matrix[i,i]}, are in i-th column of \code{Condition}. 
#' For C-vine: data corresponding to the conditioning variable whose index is in \code{RVM$Matrix[i,i]}, are in [(d+1)-i]-th column 
#' of \code{Condition}. See examples.
#' @param N Number of data to be simulated. By default N is taken from \code{Condition}, which is a \code{N x Nx} matrix.
#' It is necessary to specify \code{N} only when \code{Condition} is not given.
#'
#' @return A \code{N x d} matrix of the simulated variables from the given C- or D-vine copula model. In the first columns there are
#' the simulated conditioned variables, and in the last columns the conditioning variables \code{Condition}.
#' For more details about the exact order of the variables in the columns see the examples. The 
#' function is built to work easily in combination with \code{\link{CDVineCondFit}}. 
#'
#' 
#'
#' @examples 
#' 
#' # Example 1: conditional sampling from a C-Vine
#' 
#' # Read data
#' data(dataset) 
#' data <- dataset$data[1:400,1:4]
#'
#' # Define the variables Y and X. X are the conditioning variables, 
#' # which have to be positioned in the last columns of the data.frame
#' colnames(data) <- c("Y1","Y2","X3","X4")
#'
#' \dontrun{
#' # Select a vine and fit the copula families, specifying that there are 2 conditioning variables
#' RVM <- CDVineCondFit(data,Nx=2,type="CVine")
#'
#' # Set the values of the conditioning variables as those used for the calibration. 
#' # Order them with respect to RVM$Matrix, considering that is a C-Vine
#' d=dim(RVM$Matrix)[1]
#' cond1 <- data[,RVM$Matrix[(d+1)-1,(d+1)-1]]
#' cond2 <- data[,RVM$Matrix[(d+1)-2,(d+1)-2]]
#' condition <- cbind(cond1,cond2)
#'
#' # Simulate the variables
#' Sim <- CDVineCondSim(RVM,condition)
#'
#' # Plot the simulated variables over the observed
#' Sim <- data.frame(Sim)
#' overplot(Sim,data)
#' 
#' 
#' 
#' # Example 2: conditional sampling from a D-Vine
#'
#' # Read data
#' data(dataset) 
#' data <- dataset$data[1:100,1:4]
#'
#' # Define the variables Y and X. X are the conditioning variables, 
#' # which have to be positioned in the last columns of the data.frame
#' colnames(data) <- c("Y1","Y2","X3","X4")
#'
#' # Select a vine and fit the copula families, specifying that there are 2 conditioning variables
#' RVM <- CDVineCondFit(data,Nx=2,type="DVine")
#' summary(RVM) #It is a D-Vine.
#'
#' # Set the values of the conditioning variables as those used for the calibration. 
#' # Order them with respect to RVM$Matrix, considering that is a D-Vine.
#' cond1 <- data[,RVM$Matrix[1,1]]
#' cond2 <- data[,RVM$Matrix[2,2]]
#' condition <- cbind(cond1,cond2)
#'
#' # Simulate the variables
#' Sim <- CDVineCondSim(RVM,condition)
#'
#' # Plot the simulated variables over the observed
#' Sim <- data.frame(Sim)
#' overplot(Sim,data)
#' 
#' 
#' 
#' # Example 3
#'
#' # Read data
#' data(dataset) 
#' data <- dataset$data[1:100,1:2]
#' colnames(data) <- c("Y1","X2")
#'
#' # Fit copula
#' require(VineCopula)
#' BiCop <- BiCopSelect(data$Y1,data$X2)
#' BiCop
#'
#' # Fix conditioning variable to low values and simulate
#' condition <- data$X2/10
#' Sim <- CDVineCondSim(BiCop,condition)
#'
#' # Plot the simulated variables over the observed
#' Sim <- data.frame(Sim)
#' overplot(Sim,data)
#' }
#' 
#' @author Emanuele Bevacqua
#'
#' @references Bevacqua, E., Maraun, D., Hobaek Haff, I., Widmann, M., and Vrac, M.: Multivariate statistical modelling of compound events via pair-copula constructions: analysis of floods in Ravenna (Italy), 
#' Hydrol. Earth Syst. Sci., 21, 2701-2723, https://doi.org/10.5194/hess-21-2701-2017, 2017.
#' \href{https://www.researchgate.net/publication/317414374_Multivariate_statistical_modelling_of_compound_events_via_pair-copula_constructions_Analysis_of_floods_in_Ravenna_Italy}{[link]} 
#' \href{https://www.hydrol-earth-syst-sci.net/21/2701/2017/hess-21-2701-2017.html}{[link]} 
#'
#' Aas, K., Czado, C., Frigessi, A. and Bakken, H.: Pair-copula constructions of multiple dependence, Insurance:
#' Mathematics and Economics, 44(2), 182-198, <doi:10.1016/j.insmatheco.2007.02.001>, 2009. 
#' \href{http://www.sciencedirect.com/science/article/pii/S0167668707000194}{[link]} 
#'
#' Ulf Schepsmeier, Jakob Stoeber, Eike Christian Brechmann, Benedikt Graeler, Thomas 
#' Nagler and Tobias Erhardt (2017). VineCopula: Statistical Inference of Vine Copulas. R 
#' package version 2.1.1. \href{https://CRAN.R-project.org/package=VineCopula}{[link]}
#' 
#' @seealso \code{\link{CDVineCondFit}}
#' @import VineCopula
#' @export CDVineCondSim
CDVineCondSim <- function(RVM,Condition,N)
{
  d <- dim(as.matrix(RVM$family))[1]
  if(d==1)
  {
    return(DVineCondSim(RVM,Condition,N))
  }
  if(is.matrix(Condition))
     {
       if(dim(Condition)[2]==d)
         {
           print("Please, provide a number of conditioning variables smaller than d")
         }
      }
  
  
  if(RVM$Matrix[d,1]!=RVM$Matrix[d,2])
  {
    return(DVineCondSim(RVM,Condition,N))
  }
  else if(RVM$Matrix[d,1]==RVM$Matrix[d,2])
  {
    return(CVineCondSim(RVM,Condition,N))
  }
}

























#' Ranking of C- and D- vines allowing for conditional simulation
#'
#' @description Provides a ranking of the C- and D- vines which allow for conditional 
#' sampling, under the condition discussed in the descriprion of \code{\link{CDVineCondFit}}.
#'
#' @param data An \code{N x d} data matrix (with uniform margins).
#' Data of the conditioning variable(s) have to occupy the last column(s) of this matrix.
#'
#' @param Nx Number of conditioning variables.
#'
#' @param treecrit Character indicating the criteria used to select the vine. All possible vines are fitted trough 
#' the function \code{\link[VineCopula]{RVineCopSelect}} of the package \code{VineCopula}. Then the vines are ranked with respect 
#' the Akaike information criterion (\code{treecrit = "AIC"}, default) or Bayesian information criterion (\code{treecrit = "BIC"}). 
#' This need the estimation and model selection for all the pairs of all the possible vines, therefore could 
#' require long time in case of large datasets, i.e. large \code{N x d}.
#' 
#' @param selectioncrit Character indicating the criterion for pair-copula selection.
#' Possible choices are \code{selectioncrit = "AIC"} (default) and \code{"BIC"}.
#'
#' @param familyset "Integer vector of pair-copula families to select from. The vector has to include at least one
#' pair-copula family that allows for positive and one that allows for negative
#' dependence. Not listed copula families might be included to better handle
#' limit cases.  If \code{familyset = NA} (default), selection among all
#' possible families is performed. If a vector of negative numbers is provided,
#' selection among all but \code{abs(familyset)} is performed. Coding of bivariate copula families: \cr
#' \code{0} = independence copula \cr
#' \code{1} = Gaussian copula \cr
#' \code{2} = Student t copula (t-copula) \cr
#' \code{3} = Clayton copula \cr
#' \code{4} = Gumbel copula \cr
#' \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr
#' \code{7} = BB1 copula \cr
#' \code{8} = BB6 copula \cr
#' \code{9} = BB7 copula \cr
#' \code{10} = BB8 copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' \code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' \code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' \code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' \code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr
#' \code{24} = rotated Gumbel copula (90 degrees) \cr
#' \code{26} = rotated Joe copula (90 degrees) \cr
#' \code{27} = rotated BB1 copula (90 degrees) \cr
#' \code{28} = rotated BB6 copula (90 degrees) \cr
#' \code{29} = rotated BB7 copula (90 degrees) \cr
#' \code{30} = rotated BB8 copula (90 degrees) \cr
#' \code{33} = rotated Clayton copula (270 degrees) \cr
#' \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr
#' \code{37} = rotated BB1 copula (270 degrees) \cr
#' \code{38} = rotated BB6 copula (270 degrees) \cr
#' \code{39} = rotated BB7 copula (270 degrees) \cr
#' \code{40} = rotated BB8 copula (270 degrees) \cr
#' \code{104} = Tawn type 1 copula \cr
#' \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
#' \code{124} = rotated Tawn type 1 copula (90 degrees) \cr
#' \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr
#' \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
#' \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees)" (VineCopula Documentation, version 2.1.1, pp. 73-74)
#'
#' 
#'
#' @param type Type of vine to be fitted: \cr
#' C-Vine: "CVine" or 1; \cr
#' D-Vine: "DVine" or 2; \cr
#' Both C and D-Vine: "CVine-DVine" or "1-2" (default).
#'
#' @param indeptest "Logical; whether a hypothesis test for the independence of
#' \code{u1} and \code{u2} is performed before bivariate copula selection
#' (default: \code{indeptest = FALSE}; see \code{BiCopIndTest}).  The
#' independence copula is chosen for a (conditional) pair if the null
#' hypothesis of independence cannot be rejected.
#'
#' @param level numeric; significance level of the independence test (default:
#' \code{level = 0.05}). 
#'
#' @param se Logical; whether standard errors are estimated (default: \code{se
#' = FALSE}). 
#'
#' @param rotations logical; if \code{TRUE}, all rotations of the families in
#' \code{familyset} are included.
#'
#' @param method indicates the estimation method: either maximum
#' likelihood estimation (\code{method = "mle"}; default) or inversion of
#' Kendall's tau (\code{method = "itau"}). For \code{method = "itau"} only
#' one parameter families and the Student t copula can be used (\code{family =
#' 1,2,3,4,5,6,13,14,16,23,24,26,33,34} or \code{36}). For the t-copula,
#' \code{par2} is found by a crude profile likelihood optimization over the
#' interval (2, 10]." (VineCopula Documentation, version 2.1.1, pp. 74-75)
#'
#' @return \describe{
#' \item{\code{table}}{A table with the ranking of the vines, with vine index \code{i},
#' values of the selected \code{treecrit} and vine \code{type} (1 for "CVine" and 2 for D-Vine).}
#'
#' \item{\code{vines}}{A list where the element \code{[[i]]} is an \code{\link[VineCopula]{RVineMatrix}} object corresponding to
#' the \code{i}-th vine in the ranking shown in \code{table}. 
#' Each \code{\link[VineCopula]{RVineMatrix}} object containes the selected families (\code{$family}) as well as sequentially
#' estimated parameters stored in \code{$par} and \code{$par2}. Details about \code{\link[VineCopula]{RVineMatrix}} objects are given in 
#' the documentation file of the \code{VineCopula} package). 
#' The fit of each model is performed via the function \code{\link[VineCopula]{RVineCopSelect}} of the package \code{VineCopula}. 
#' "The object is augmented by the following information about the fit: 
#' 
#' \describe{
#' \item{\code{se, se2}}{standard errors for the parameter estimates  (if \code{se = TRUE}; note that these are only 
#' approximate since they do not account for the sequential nature of the estimation}
#'
#' \item{\code{nobs}}{number of observations}
#'
#' \item{\code{logLik, pair.logLik}}{log likelihood (overall and pairwise)}
#'
#' \item{\code{AIC, pair.AIC}}{Aikaike's Informaton Criterion (overall and pairwise)}
#'
#' \item{\code{BIC, pair.BIC}}{Bayesian's Informaton Criterion (overall and pairwise)}
#'
#' \item{\code{emptau}}{matrix of empirical values of Kendall's tau}
#'
#' \item{\code{p.value.indeptest}}{matrix of p-values of the independence test.}}
#' }}
#'
#' @note For a comprehensive summary of the vine copula model, use
#' \code{summary(object)}; to see all its contents, use \code{str(object)}." (VineCopula Documentation, version 2.1.1, pp. 103)
#'
#' @examples
#' 
#' # Read data
#' data(dataset)
#' data <- dataset$data[1:100,1:5]
#'
#' # Define the variables Y and X. X are the conditioning variables, 
#' # which have to be positioned in the last columns of the data.frame
#' colnames(data) <- c("Y1","Y2","X3","X4","X5")
#'
#' # Rank the possible D-Vines according to the AIC
#' \dontrun{
#' Ranking <- CDVineCondRank(data,Nx=3,"AIC",type="DVine")
#' Ranking$table
#' #    tree       AIC type
#' # 1     1 -292.8720    2
#' # 2     2 -290.2941    2
#' # 3     3 -288.5719    2
#' # 4     4 -288.2496    2
#' # 5     5 -287.8006    2
#' # 6     6 -285.8503    2
#' # 7     7 -282.2867    2
#' # 8     8 -278.9371    2
#' # 9     9 -275.8339    2
#' # 10   10 -272.9459    2
#' # 11   11 -271.1526    2
#' # 12   12 -270.5269    2
#' 
#' Ranking$vines[[1]]$AIC
#' # [1] -292.8720
#' summary(Ranking$vines[[1]])
#' }
#' 
#' @author Emanuele Bevacqua
#' 
#' @references Bevacqua, E., Maraun, D., Hobaek Haff, I., Widmann, M., and Vrac, M.: Multivariate statistical modelling of compound events via pair-copula constructions: analysis of floods in Ravenna (Italy), 
#' Hydrol. Earth Syst. Sci., 21, 2701-2723, https://doi.org/10.5194/hess-21-2701-2017, 2017.
#' \href{https://www.researchgate.net/publication/317414374_Multivariate_statistical_modelling_of_compound_events_via_pair-copula_constructions_Analysis_of_floods_in_Ravenna_Italy}{[link]} 
#' \href{https://www.hydrol-earth-syst-sci.net/21/2701/2017/hess-21-2701-2017.html}{[link]} 
#' 
#' Aas, K., Czado, C., Frigessi, A. and Bakken, H.: Pair-copula constructions of multiple dependence, Insurance:
#' Mathematics and Economics, 44(2), 182-198, <doi:10.1016/j.insmatheco.2007.02.001>, 2009. 
#' \href{http://www.sciencedirect.com/science/article/pii/S0167668707000194}{[link]} 
#' 
#' Ulf Schepsmeier, Jakob Stoeber, Eike Christian Brechmann, Benedikt Graeler, Thomas 
#' Nagler and Tobias Erhardt (2017). VineCopula: Statistical Inference of Vine Copulas. R 
#' package version 2.1.1. \href{https://CRAN.R-project.org/package=VineCopula}{[link]}
#' 
#' @seealso \code{\link{CDVineCondFit}}
#' @import VineCopula combinat
#' @export CDVineCondRank
CDVineCondRank <- function(data,Nx,treecrit="AIC",selectioncrit="AIC",familyset = NA,type="CVine-DVine",indeptest = FALSE, level = 0.05, se = FALSE,
                           rotations = TRUE, method = "mle"){
  if(type==1)
  {
    type="CVine"
  }
  if(type==2)
  {
    type="DVine"
  }
  if(type=="1-2")
  {
    type="CVine-DVine"
  }

  if(type=="CVine-DVine" | type=="CVine")
  {
    RankingC <- RankCVineCond(data,Nx,treecrit,selectioncrit,familyset,indeptest = indeptest, level = level, se = se,
                              rotations = rotations, method = method)
  }
  if(type=="CVine-DVine" | type=="DVine")
  {
    RankingD <- RankDVineCond(data,Nx,treecrit,selectioncrit,familyset,indeptest = indeptest, level = level, se = se,
                              rotations = rotations, method = method)
  }
  
    if(type=="CVine-DVine")
    {
      RankingCD=list()
      RankingCD[[1]]=data.frame(tree=seq(1,(length(RankingD[[1]]$tree)+length(RankingC[[1]]$tree))),
                                AICBIC=rep(-9999,(length(RankingD[[1]]$tree)+length(RankingC[[1]]$tree))),
                                type=rep(-9999,(length(RankingD[[1]]$tree)+length(RankingC[[1]]$tree))))
      RankingCD[[2]]=list()
      indexC=1
      indexD=1
      for(i in 1:(length(RankingD[[1]]$tree)+length(RankingC[[1]]$tree)))
      {
        if(indexC<=length(RankingC[[1]]$tree) & indexD<=length(RankingD[[1]]$tree))
        {
          if(RankingC[[1]][indexC,2]<RankingD[[1]][indexD,2])
          {
            RankingCD[[1]]$AICBIC[i]=RankingC[[1]][indexC,2]
            RankingCD[[1]]$type[i]=1
            RankingCD[[2]][[i]]=RankingC[[2]][[indexC]]
            indexC=indexC+1
          }
          else
          {
            RankingCD[[1]]$AICBIC[i]=RankingD[[1]][indexD,2]
            RankingCD[[1]]$type[i]=2
            RankingCD[[2]][[i]]=RankingD[[2]][[indexD]]
            indexD=indexD+1
          }
        }
        else if(indexC>length(RankingC[[1]]$tree))
        {
          RankingCD[[1]]$AICBIC[i]=RankingD[[1]][indexD,2]
          RankingCD[[1]]$type[i]=2
          RankingCD[[2]][[i]]=RankingD[[2]][[indexD]]
          indexD=indexD+1
        }
        else if(indexD>length(RankingD[[1]]$tree))
        {
          RankingCD[[1]]$AICBIC[i]=RankingC[[1]][indexC,2]
          RankingCD[[1]]$type[i]=1
          RankingCD[[2]][[i]]=RankingC[[2]][[indexC]]
          indexC=indexC+1
        }
      }

      if(treecrit=="AIC")
      {
        colnames(RankingCD[[1]])=c("tree","AIC","type")
      }
      else if(treecrit=="BIC")
      {
        colnames(RankingCD[[1]])=c("tree","BIC","type")
      }
      Rank=list()
      Rank$table=RankingCD[[1]]
      Rank$vines=RankingCD[[2]]
      return(Rank)
    }
    else if(type=="CVine")
    {
      RankingC[[1]]$type=1
      Rank=list()
      Rank$table=RankingC[[1]]
      Rank$vines=RankingC[[2]]
      return(Rank)
    }
    else if(type=="DVine")
    {
      RankingD[[1]]$type=2
      Rank=list()
      Rank$table=RankingD[[1]]
      Rank$vines=RankingD[[2]]
      return(Rank)
    }
}


























#' List of the possible C- and D- vines allowing for conditional simulation
#'
#' @description Provides a list of the C- and D- vines which allow for conditional 
#' sampling, under the condition discussed in the descriprion of \code{\link{CDVineCondFit}}.
#'
#' @param data An \code{N x d} data matrix.
#' Data of the conditioning variable(s) have to occupy the last column(s) of this matrix.
#'
#' @param Nx Number of conditioning variables.
#'
#' @param type Type of vine to be considered: \cr
#' C-Vine: "CVine" or 1; \cr
#' D-Vine: "DVine" or 2; \cr
#' Both C and D-Vine: "CVine-DVine" or "1-2" (default).
#'
#' @return Listes of matrices describing C- (\code{$CVine}) and D- (\code{$DVine}) Vines. 
#' Each matrix corresponds to a vine, according to the same notation used for \code{\link[VineCopula]{RVineMatrix}} 
#' objects (for further details about \code{\link[VineCopula]{RVineMatrix}} objects see the documentation file of the \code{VineCopula} package). 
#' The index \code{i} in the matrix corresponds to the variable in the i-th column of \code{data}. 
#' 
#' @examples
#' 
#' # Read data
#' data(dataset)
#' data <- dataset$data[1:100,1:5]
#'
#' # Define the variables Y and X. X are the conditioning variables, 
#' # which have to be positioned in the last columns of the data.frame
#' colnames(data) <- c("Y1","Y2","X3","X4","X5")
#'
#' # List possible D-Vines:
#' ListVines <- CDVineCondListMatrices(data,Nx=3,"DVine")
#' ListVines$DVine
#' 
#' # List possible C-Vines:
#' ListVines <- CDVineCondListMatrices(data,Nx=3,"CVine")
#' ListVines$CVine
#' 
#' # List possible C- and D-Vines:
#' ListVines <- CDVineCondListMatrices(data,Nx=3,"CVine-DVine")
#' ListVines
#' 
#' @author Emanuele Bevacqua
#' 
#' @references Bevacqua, E., Maraun, D., Hobaek Haff, I., Widmann, M., and Vrac, M.: Multivariate statistical modelling of compound events via pair-copula constructions: analysis of floods in Ravenna (Italy), 
#' Hydrol. Earth Syst. Sci., 21, 2701-2723, https://doi.org/10.5194/hess-21-2701-2017, 2017.
#' \href{https://www.researchgate.net/publication/317414374_Multivariate_statistical_modelling_of_compound_events_via_pair-copula_constructions_Analysis_of_floods_in_Ravenna_Italy}{[link]} 
#' \href{https://www.hydrol-earth-syst-sci.net/21/2701/2017/hess-21-2701-2017.html}{[link]} 
#'
#' Aas, K., Czado, C., Frigessi, A. and Bakken, H.: Pair-copula constructions of multiple dependence, Insurance:
#' Mathematics and Economics, 44(2), 182-198, <doi:10.1016/j.insmatheco.2007.02.001>, 2009. 
#' \href{http://www.sciencedirect.com/science/article/pii/S0167668707000194}{[link]} 
#' 
#' @seealso \code{\link{CDVineCondFit}}
#' @import VineCopula combinat
#' @export CDVineCondListMatrices
CDVineCondListMatrices <- function(data,Nx,type="CVine-DVine"){
  if(type==1)
  {
    type="CVine"
  }
  if(type==2)
  {
    type="DVine"
  }
  if(type=="1-2")
  {
    type="CVine-DVine"
  }
  
  
  if(type=="CVine-DVine" | type=="CVine")
  {
    MatricesC <- PossibleCVineMatrixCond(data,Nx)
    if(type=="CVine")
    {
      a=list()
      a$CVine=MatricesC
      return(a)
    }
  }
  
  if(type=="CVine-DVine" | type=="DVine")
  {
    MatricesD <- PossibleDVineMatrixCond(data,Nx)
    if(type=="DVine")
    {
      a=list()
      a$DVine=MatricesD
      return(a)
    }
  }
  if(type=="CVine-DVine")
  {
    a=list()
    a$CVine=MatricesC
    a$DVine=MatricesD
    return(a)
  }
}





















#' Selection of a C- or D- vine copula model for conditional sampling
#'
#' @description This function fits either a C- or a D- vine model to a d-dimensional dataset of uniform variables. 
#' The fit of the pair-copula families is performed sequentially through the function \code{\link[VineCopula]{RVineCopSelect}} of 
#' the package \code{VineCopula}. The vine structure is selected among a group of C- and a D- vines which satisfy the requirement 
#' discussed in \emph{Bevacqua et al. (2017)}. This group is composed by all C- and D- vines from which the conditioning variables 
#' would be sampled as first when following the algorithms from \emph{Aas et al. (2009)}. Alternatively, if the 
#' vine matrix describing the vine structure is given to the function, the fit of the pair-copulas is directly performed skipping the vine structure 
#' selection procedure.
#' 
#' @param data An \code{N x d} data matrix (with uniform margins).
#' The data of the conditioning variable(s) have to occupy the last column(s) of this matrix.
#'
#' @param Nx Number of conditioning variables.
#'
#' @param treecrit Character indicating the criteria used to select the vine. All possible vines are fitted trough 
#' the function \code{\link[VineCopula]{RVineCopSelect}} of the package \code{VineCopula}. Then the vines are ranked with respect 
#' the Akaike information criterion (\code{treecrit = "AIC"}, default) or Bayesian information criterion (\code{treecrit = "BIC"}). 
#' This need the estimation and model selection for all the pairs of all the possible vines, therefore could 
#' require long time in case of large datasets, i.e. large \code{N x d}.
#' 
#' @param type Type of vine to be fitted: \cr
#' C-Vine: "CVine" or 1; \cr
#' D-Vine: "DVine" or 2; \cr
#' Both C and D-Vine: "CVine-DVine" or "1-2" (default).
#' 
#' @param selectioncrit Character indicating the criterion for pair-copula selection.
#' Possible choices are \code{"AIC"} (default) and \code{"BIC"}.
#'
#' @param familyset "Integer vector of pair-copula families to select from. The vector has to include at least one
#' pair-copula family that allows for positive and one that allows for negative
#' dependence. Not listed copula families might be included to better handle
#' limit cases.  If \code{familyset = NA} (default), selection among all
#' possible families is performed. If a vector of negative numbers is provided,
#' selection among all but \code{abs(familyset)} is performed. Coding of bivariate copula families: \cr
#' \code{0} = independence copula \cr
#' \code{1} = Gaussian copula \cr
#' \code{2} = Student t copula (t-copula) \cr
#' \code{3} = Clayton copula \cr
#' \code{4} = Gumbel copula \cr
#' \code{5} = Frank copula \cr
#' \code{6} = Joe copula \cr
#' \code{7} = BB1 copula \cr
#' \code{8} = BB6 copula \cr
#' \code{9} = BB7 copula \cr
#' \code{10} = BB8 copula \cr
#' \code{13} = rotated Clayton copula (180 degrees; ``survival Clayton'') \cr
#' \code{14} = rotated Gumbel copula (180 degrees; ``survival Gumbel'') \cr
#' \code{16} = rotated Joe copula (180 degrees; ``survival Joe'') \cr
#' \code{17} = rotated BB1 copula (180 degrees; ``survival BB1'')\cr
#' \code{18} = rotated BB6 copula (180 degrees; ``survival BB6'')\cr
#' \code{19} = rotated BB7 copula (180 degrees; ``survival BB7'')\cr
#' \code{20} = rotated BB8 copula (180 degrees; ``survival BB8'')\cr
#' \code{23} = rotated Clayton copula (90 degrees) \cr
#' \code{24} = rotated Gumbel copula (90 degrees) \cr
#' \code{26} = rotated Joe copula (90 degrees) \cr
#' \code{27} = rotated BB1 copula (90 degrees) \cr
#' \code{28} = rotated BB6 copula (90 degrees) \cr
#' \code{29} = rotated BB7 copula (90 degrees) \cr
#' \code{30} = rotated BB8 copula (90 degrees) \cr
#' \code{33} = rotated Clayton copula (270 degrees) \cr
#' \code{34} = rotated Gumbel copula (270 degrees) \cr
#' \code{36} = rotated Joe copula (270 degrees) \cr
#' \code{37} = rotated BB1 copula (270 degrees) \cr
#' \code{38} = rotated BB6 copula (270 degrees) \cr
#' \code{39} = rotated BB7 copula (270 degrees) \cr
#' \code{40} = rotated BB8 copula (270 degrees) \cr
#' \code{104} = Tawn type 1 copula \cr
#' \code{114} = rotated Tawn type 1 copula (180 degrees) \cr
#' \code{124} = rotated Tawn type 1 copula (90 degrees) \cr
#' \code{134} = rotated Tawn type 1 copula (270 degrees) \cr
#' \code{204} = Tawn type 2 copula \cr
#' \code{214} = rotated Tawn type 2 copula (180 degrees) \cr
#' \code{224} = rotated Tawn type 2 copula (90 degrees) \cr
#' \code{234} = rotated Tawn type 2 copula (270 degrees) 
#'
#'
#' @param indeptest Logical; whether a hypothesis test for the independence of
#' \code{u1} and \code{u2} is performed before bivariate copula selection
#' (default: \code{indeptest = FALSE}; see \code{BiCopIndTest}).  The
#' independence copula is chosen for a (conditional) pair if the null
#' hypothesis of independence cannot be rejected. 
#'
#' @param level numeric; significance level of the independence test (default:
#' \code{level = 0.05}). 
#'
#' @param se Logical; whether standard errors are estimated (default: \code{se
#' = FALSE}). 
#'
#' @param rotations logical; if \code{TRUE}, all rotations of the families in
#' \code{familyset} are included. 
#'
#' @param method indicates the estimation method: either maximum
#' likelihood estimation (\code{method = "mle"}; default) or inversion of
#' Kendall's tau (\code{method = "itau"}). For \code{method = "itau"} only
#' one parameter families and the Student t copula can be used (\code{family =
#' 1,2,3,4,5,6,13,14,16,23,24,26,33,34} or \code{36}). For the t-copula,
#' \code{par2} is found by a crude profile likelihood optimization over the
#' interval (2, 10]." (VineCopula Documentation, version 2.1.1, pp. 73-75)
#' 
#' @param Matrix \code{d x d} matrix that defines the vine structure. 
#' If \code{Matrix} is not given, the routine finds the best vine structure according to \code{selectioncrit}. 
#' If \code{Matrix} is given, the fit is performed only if the structure respects the necessary conditions 
#' for the conditional sampling (if the conditions are not respected, an error message is returned).
# 
#' @return An \code{\link[VineCopula]{RVineMatrix}} object describing the selected copula model 
#' (for further details about \code{\link[VineCopula]{RVineMatrix}} objects see the documentation file of the \code{VineCopula} package). 
#' The selected families are stored 
#' in \code{$family}, and the sequentially estimated parameters in \code{$par} and \code{$par2}. 
#' The fit of the model is performed via the function \code{\link[VineCopula]{RVineCopSelect}} of the package \code{VineCopula}.
#'  
#' "The object \code{\link[VineCopula]{RVineMatrix}} includes the following information about the fit:
#' \describe{
#' \item{\code{se, se2}}{standard errors for the parameter estimates  (if \code{se = TRUE}; note that these are only 
#' approximate since they do not account for the sequential nature of the estimation,}  
#'
#' \item{\code{nobs}}{number of observations,} 
#'
#' \item{\code{logLik, pair.logLik}}{log likelihood (overall and pairwise)} 
#'
#' \item{\code{AIC, pair.AIC}}{Aikaike's Informaton Criterion (overall and pairwise),} 
#'
#' \item{\code{BIC, pair.BIC}}{Bayesian's Informaton Criterion (overall and pairwise),} 
#'
#' \item{\code{emptau}}{matrix of empirical values of Kendall's tau,} 
#'
#' \item{\code{p.value.indeptest}}{matrix of p-values of the independence test.}
#' }
#'
#' @note For a comprehensive summary of the vine copula model, use
#' \code{summary(object)}; to see all its contents, use \code{str(object)}". (VineCopula Documentation, version 2.1.1, pp. 103)
#'
#' @examples
#' 
#' # Example 1
#' 
#' # Read data 
#' data(dataset) 
#' data <- dataset$data[1:100,1:5]
#' 
#' # Define the variables Y and X. X are the conditioning variables, 
#' # which have to be positioned in the last columns of the data.frame
#' colnames(data) <- c("Y1","Y2","X3","X4","X5")
#'
#' \dontrun{
#' # Select and fit a C- vine copula model, requiring that the 
#' RVM <- CDVineCondFit(data,Nx=3,treecrit="BIC",type="CVine",selectioncrit="AIC")
#' summary(RVM)
#' RVM$Matrix
#' }
#' 
#' 
#' 
#' # Example 2
#' 
#' # Read data 
#' data(dataset) 
#' data <- dataset$data[1:80,1:5]
#' 
#' # Define the variables Y and X. X are the conditioning variables, 
#' # which have to be positioned in the last columns of the data.frame
#' colnames(data) <- c("Y1","Y2","X3","X4","X5")
#'
#' # Define a VineMatrix which can be used for conditional sampling
#' ListVines <- CDVineCondListMatrices(data,Nx=3)
#' Matrix=ListVines$DVine[[1]]
#' Matrix
#' 
#' \dontrun{
#' # Fit copula families for the defined vine:
#' RVM <- CDVineCondFit(data,Nx=3,Matrix=Matrix)
#' summary(RVM)
#' RVM$Matrix
#' RVM$family
#' 
#' # check
#' identical(RVM$Matrix,Matrix)
#' 
#' # Fit copula families for the defined vine, given a group of families to select from:
#' RVM <- CDVineCondFit(data,Nx=3,Matrix=Matrix,familyset=c(1,2,3,14))
#' summary(RVM)
#' RVM$Matrix
#' RVM$family
#' 
#' # Try to fit copula families for a vine which is not among those 
#' # that allow for conditional sampling:
#' Matrix
#' Matrix[which(Matrix==4)]=40
#' Matrix[which(Matrix==2)]=20
#' Matrix[which(Matrix==40)]=2
#' Matrix[which(Matrix==20)]=4
#' Matrix
#' RVM <- CDVineCondFit(data,Nx=3,Matrix=Matrix)
#' RVM
#' }
#'
#' @author Emanuele Bevacqua
#' 
#' @references Bevacqua, E., Maraun, D., Hobaek Haff, I., Widmann, M., and Vrac, M.: Multivariate statistical modelling of compound events via pair-copula constructions: analysis of floods in Ravenna (Italy), 
#' Hydrol. Earth Syst. Sci., 21, 2701-2723, https://doi.org/10.5194/hess-21-2701-2017, 2017.
#' \href{https://www.researchgate.net/publication/317414374_Multivariate_statistical_modelling_of_compound_events_via_pair-copula_constructions_Analysis_of_floods_in_Ravenna_Italy}{[link]} 
#' \href{https://www.hydrol-earth-syst-sci.net/21/2701/2017/hess-21-2701-2017.html}{[link]} 
#' 
#' Aas, K., Czado, C., Frigessi, A. and Bakken, H.: Pair-copula constructions of multiple dependence, Insurance: 
#' Mathematics and Economics, 44(2), 182-198, <doi:10.1016/j.insmatheco.2007.02.001>, 2009.
#' \href{http://www.sciencedirect.com/science/article/pii/S0167668707000194}{[link]} 
#'
#' Ulf Schepsmeier, Jakob Stoeber, Eike Christian Brechmann, Benedikt Graeler, Thomas 
#' Nagler and Tobias Erhardt (2017). VineCopula: Statistical Inference of Vine Copulas. R 
#' package version 2.1.1. \href{https://CRAN.R-project.org/package=VineCopula}{[link]}
#' 
#' @seealso \code{\link{CDVineCondSim}}, \code{\link{CDVineCondRank}}
#' @export CDVineCondFit
CDVineCondFit <- function(data,Nx,treecrit="AIC",type="CVine-DVine",selectioncrit="AIC",familyset=NA,indeptest = FALSE, level = 0.05, se = FALSE,
                          rotations = TRUE, method = "mle",Matrix=FALSE)
{
  
  if(!is.numeric(familyset)[1])
  {
    familyset=c(1,2,3,4,5,6,7,8,9,10,13,14,16,17,18,19,20,23,24,26,27,28,29,30,33,34,36,37,38,39,
                40,104,114,124,134,
                204,
                214,
                224,
                234
    )
  }
  
  if(!is.null(dim(Matrix)))
  {
    #quality check
    ListVines <- CDVineCondListMatrices(data,Nx,"CVine-DVine")
    quality=FALSE
    for(i in 1:2)
      {
       for(j in 1:length(ListVines[[i]]))
          {
             if(identical(ListVines[[i]][[j]],Matrix))
              {
                quality=TRUE
                break
              }
          }
      }
      
    if(quality==TRUE)
      {
        vineFIT <- RVineCopSelect(data,familyset,Matrix,selectioncrit,indeptest = indeptest, level = level, se = se,
                              rotations = rotations, method = method)
        return(vineFIT)
      }
    else
     {
      print("The given Matrix cannot be used for the conditional simulation. To fit a generic Vine, please use the function RvineCopSelect from the package VineCopula")
      return(NULL)
    }
  }
  
  
  if(type==1)
  {
    type="CVine"
  }
  if(type==2)
  {
    type="DVine"
  }
  if(type=="1-2")
  {
    type="CVine-DVine"
  }

  
  
  if(type=="CVine-DVine" | type=="CVine")
  {
    MatrixC <- FirstCVineCond(data,Nx=Nx,treecrit=treecrit,selectioncrit=selectioncrit,familyset=familyset,indeptest = indeptest, level = level, se = se,
                              rotations = rotations, method = method)
  }
  if(type=="CVine-DVine" | type=="DVine")
  {
    MatrixD <- FirstDVineCond(data,Nx=Nx,treecrit=treecrit,selectioncrit=selectioncrit,familyset=familyset,indeptest = indeptest, level = level, se = se,
                              rotations = rotations, method = method)
  }

  
    if(type=="CVine-DVine")
    {
      if(MatrixC[[2]]<MatrixD[[2]])
      {
        vineFIT <- RVineCopSelect(data,familyset,MatrixC[[1]],selectioncrit,indeptest = indeptest, level = level, se = se,
                                  rotations = rotations, method = method)
      }
      else
      {
        vineFIT <- RVineCopSelect(data,familyset,MatrixD[[1]],selectioncrit,indeptest = indeptest, level = level, se = se,
                                  rotations = rotations, method = method)
      }
      return(vineFIT)
    }
    else if(type=="CVine")
    {
      vineFIT <- RVineCopSelect(data,familyset,MatrixC[[1]],selectioncrit,indeptest = indeptest, level = level, se = se,
                                rotations = rotations, method = method)
      return(vineFIT)
    }
    else if(type=="DVine")
    {
      vineFIT <- RVineCopSelect(data,familyset,MatrixD[[1]],selectioncrit,indeptest = indeptest, level = level, se = se,
                                rotations = rotations, method = method)
      return(vineFIT)
    }
}


























#' overplot
#' @description This function overlays the scatterplot matrices of two multivariate datsets. 
#' Moreover, it shows the dependencies among all the pairs for both datsets.
#'
#' @param data1,data2 Two \code{N x d} matrices of data to be plotted.
#' @param col1,col2 Colors used for \code{data1} and \code{data2} during the plot. Default is \code{col1="black"} and \code{col2="grey"}.
#' @param pch1,pch2 Paramter to specify the symbols to use when plotting points of \code{data1} and \code{data2}.Default is \code{pch1=1} and \code{pch2=1}.
#' @param xlim,ylim Two bidimensional vectors indicating the limits of x and y axes for all the scatterplots. 
#' If not given, they are authomatically 
#' computed for each of the scatterplots.
#' @param labels A character vector with the variable names to be printed (if not given, the names of \code{data1} 
#' variables are printed).
#' @param method Character indicating the dependence types to be computed between the pairs. Possibilites: "kendall",
#' "spearman" and "pearson" (default)
#' @param cex.cor Number: character dimension of the printed dependencies. Default \code{cex.cor=1}.
#' @param cex.labels Number: character dimension of the printed variable names. Default \code{cex.labels=1}.
#' @param cor.signif Number: number of significant numbers of the printed dependencies. Default \code{cor.signif=2}.
#' @param cex.axis Number: dimension of the axis numeric values. Default cex.axis=1.
#'
#' @return A matrix of overlaying scatterplots of the multivariate datsets \code{data1} and \code{data2}, with 
#' the dependencies of the pairs.
#' 
#' @examples
#' 
#' # Example 1
#' 
#' # Read and prepare the data for the plot
#' data(dataset) 
#' data1 <- dataset$data[1:300,]
#' data2 <- dataset$data[301:600,]
#' overplot(data1,data2,xlim=c(0,1),ylim=c(0,1),method="kendall")
#' 
#' 
#' 
#' \dontrun{
#' # Example 2
#' 
#' # Read and prepare the data for the plot
#' data(dataset) 
#' data <- dataset$data[1:200,1:5]
#' colnames(data) <- c("Y1","Y2","X3","X4","X5")
#' 
#' # Fit copula families for the defined vine:
#' ListVines <- CDVineCondListMatrices(data,Nx=3)
#' Matrix=ListVines$CVine[[1]]
#' RVM <- CDVineCondFit(data,Nx=3,Matrix=Matrix)
#' 
#' # Simulate data:
#' d=dim(RVM$Matrix)[1]
#' cond1 <- data[,RVM$Matrix[(d+1)-1,(d+1)-1]]
#' cond2 <- data[,RVM$Matrix[(d+1)-2,(d+1)-2]]
#' cond3 <- data[,RVM$Matrix[(d+1)-3,(d+1)-3]]
#' condition <- cbind(cond1,cond2,cond3)
#' Sim <- CDVineCondSim(RVM,condition)
#'
#' # Plot the simulated variables Sim over the observed
#' Sim <- data.frame(Sim)
#' overplot(data[,1:2],Sim[,1:2],xlim=c(0,1),ylim=c(0,1),method="spearman")
#' overplot(data,Sim,xlim=c(0,1),ylim=c(0,1),method="spearman")
#' }
#'
#' @author Emanuele Bevacqua
#' @importFrom stats cor
#' @importFrom graphics axis box par plot points text
#' @export overplot
overplot <- function(data1,data2,col1="black",col2="grey",xlim=NA,ylim=NA,labels=NA,method="pearson",cex.cor=1,cex.labels=1,
                     cor.signif=2,cex.axis=1,pch1=1,pch2=1){
  par(mfrow=c((length(data1[1,])-1+1),(length(data1[1,]))-1+1), oma = c(0.1, 0.1, 0.2, 0.2),
      oma = c(0,0,0,0) + 2.1,
      mar = c(0,0,1,1) + 0.1)
  
  if(!is.character(labels[1]))
  {
    labels <- colnames(data1)
  }
  for(i in 1:(length(data1[1,])-1))
  {
    plot(0,0,col="white",xaxt='n',yaxt='n')
    text(0,0,paste(labels[i]),cex=cex.labels)
    for(j in (i+1):(length(data1[1,])))
    {
      if(!is.numeric(xlim)[1])
      {
        vectJ=c(data1[,j],data2[,j])
        xlim=c(min(vectJ),max(vectJ))
      }
      if(!is.numeric(ylim)[1])
      {
        vectI=c(data1[,i],data2[,i])
        ylim=c(min(vectI),max(vectI))
      }

      plot(data1[,j],data1[,i],
           col=col1,
           xlim=xlim,ylim=ylim,
           axes=FALSE,
           xlab="",ylab="",pch=pch1)
      points(data2[,j],data2[,i],col=col2,pch=pch2)

      if(i==1 & j !=(length(data1[1,])))
      {
        axis(3,cex.axis=cex.axis)
      }
      else if(i==1 & j==(length(data1[1,])))
      {
        axis(3,cex.axis=cex.axis)
        axis(4,cex.axis=cex.axis)
      }
      else if(i!=1 & j ==(length(data1[1,])))
      {
        axis(4,cex.axis=cex.axis)
      }
      box()
    }
    if(i !=(length(data1[1,])-1))
    {
      for(k in 1:i)
      {
        plot(0,0,col="white",xaxt='n',yaxt='n')
        #plot(1, type="n", axes=F, xlab="", ylab="")
        text(0,0.2,signif(cor(data1[,k],data1[,i+1],method=method),cor.signif),col=col1,cex=cex.cor)
        text(0,-0.2,signif(cor(data2[,k],data2[,i+1],method=method),cor.signif),col=col2,cex=cex.cor)
      }
    }
  }

  for(k in 1:(length(data1[1,])-1))
  {
    plot(0,0,col="white",xaxt='n',yaxt='n')
    #plot(1, type="n", axes=F, xlab="", ylab="")
    text(0,0.2,signif(cor(data1[,k],data1[,length(data1[1,])],method=method),cor.signif),col=col1,cex=cex.cor)
    text(0,-0.2,signif(cor(data2[,k],data2[,length(data2[1,])],method=method),cor.signif),col=col2,cex=cex.cor)
  }
  plot(0,0,col="white",xaxt='n',yaxt='n')
  text(0,0,paste(labels[i+1]),cex=cex.labels)
  par(mfrow=c(1,1))
}


























































CVineCondSim <- function(RVM,Condition,N)
{
  RVM$family <- as.matrix(RVM$family)
  
  RVMFamilySim=RVM$family*0-9999
  RVMFamilySim
  #asymmetric families changed with simmetric for the simulation
  for(family in c(23,24,26,27,28,29,30))
  {
    RVMFamilySim[which(RVM$family==family)]=family+10
  }
  for(family in c(104,114,124,134))
  {
    RVMFamilySim[which(RVM$family==family)]=family+100
  }
  for(family in (10+c(23,24,26,27,28,29,30)))
  {
    RVMFamilySim[which(RVM$family==family)]=family-10
  }
  for(family in (100+c(104,114,124,134)))
  {
    RVMFamilySim[which(RVM$family==family)]=family-100
  }
  RVMFamilySim[which(RVMFamilySim==-9999)]=RVM$family[which(RVMFamilySim==-9999)]
  
  d <- dim(RVM$family)[1]
  if(d==1){d <- 2}

  w <- list()
  x <- list()
  v <- list()
  v[[1]] <- list()

  #Condition is a matrix of conditional variables
  if(missing(Condition))
  {
    Ncond=0
  }
  else if(is.vector(Condition))
  {
    cond <- list()
    cond[[1]] <- Condition
    Ncond <- 1
    N <- length(cond[[1]])
  }
  else
  {
    cond <- list()
    for(i in 1:length(Condition[1,]))
    {
      cond[[i]] <- Condition[,i]
    }
    Ncond <- length(cond)
    N <- length(cond[[1]])
  }



  if(Ncond==0)
  {
    for(t in 1:d)
    {
      w[[t]] <- BiCopSim(N,0,0)[,1]#runif(N,0,1)
    }
  }
  else if(Ncond>=1)
  {
    for(t in 1:Ncond)
    {
      w[[t]] <- cond[[t]]
    }
    for(t in ((Ncond+1):d))
    {
      w[[t]] <- BiCopSim(N,0,0)[,1]#runif(N,0,1)
    }
  }
  x[[1]] <-  w[[1]]
  v[[1]][[1]] <- x[[1]]



  for(i in 2:d)
  {
    v[[i]] <- list()
    v[[i]][[1]] <- w[[i]]

    if(i>Ncond)
    {
      for(k in (i-1):1)
      {
        v[[i]][[1]] <- BiCopHinv2(v[[i]][[1]],
                                  v[[k]][[k]],
                                  BiCop(RVMFamilySim[d-k+1,d-k+1-(i-k)],###
                                        RVM$par[d-k+1,d-k+1-(i-k)],
                                        RVM$par2[d-k+1,d-k+1-(i-k)]))#$hinv2
      }
    }

    x[[i]] <- v[[i]][[1]]
    if(i == d)
    {break}#goes out of the for(i in 2:d)
    else
    {
      for(j in 1:(i-1))
      {
        v[[i]][[j+1]] <- BiCopHfunc2(
          v[[i]][[j]],
          v[[j]][[j]],
          RVMFamilySim[d-j+1,d-j+1-(i-j)],###
          RVM$par[d-j+1,d-j+1-(i-j)],
          RVM$par2[d-j+1,d-j+1-(i-j)])#$hfunc2
      }
    }
  }

  xORDERED <- x
  matrix0 <- rep(0,d)
  for(k in 1:d)
  {
    matrix0[k] <- RVM$Matrix[k,k]
  }
  matrix0 <- rev(matrix0)
  for(k in 1:d)
  {
    xORDERED[[k]] <-x[[which(matrix0==k)]]
  }
  x <- xORDERED

  xxx <- matrix(0,nrow=N,ncol=d)
  for(k in 1:d)
  {
    xxx[,k] <- x[[k]]
  }
  return(xxx)
}









DVineCondSim <- function(RVM,Condition,N)
{
  RVM$family <- as.matrix(RVM$family)
  
  RVMFamilySim=RVM$family*0-9999
  RVMFamilySim
  #asymmetric families changed with simmetric for the simulation in BiCopHfunc2
  for(family in c(23,24,26,27,28,29,30))
  {
    RVMFamilySim[which(RVM$family==family)]=family+10
  }
  for(family in c(104,114,124,134))
  {
    RVMFamilySim[which(RVM$family==family)]=family+100
  }
  for(family in (10+c(23,24,26,27,28,29,30)))
  {
    RVMFamilySim[which(RVM$family==family)]=family-10
  }
  for(family in (100+c(104,114,124,134)))
  {
    RVMFamilySim[which(RVM$family==family)]=family-100
  }
  RVMFamilySim[which(RVMFamilySim==-9999)]=RVM$family[which(RVMFamilySim==-9999)]
  
  d <- dim(RVM$family)[1]
  if(d==1){d <- 2}

  w <- list()
  x <- list()
  v <- list()
  v[[1]] <- list()
  v[[2]] <- list()

  #Condition is a matrix of conditional variables
  if(missing(Condition))
  {
    Ncond=0
  }
  else if(is.vector(Condition))
  {
    cond <- list()
    cond[[1]] <- Condition
    Ncond <- 1
    N <- length(cond[[1]])
  }
  else
  {
    cond <- list()
    for(i in 1:length(Condition[1,]))
    {
      cond[[i]] <- Condition[,i]
    }
    Ncond <- length(cond)
    N <- length(cond[[1]])
  }



  if(Ncond==0)
  {
    for(t in 1:d)
    {
      w[[t]] <- BiCopSim(N,0,0)[,1]
    }
    x[[1]] <-  w[[1]]
    v[[1]][[1]] <- x[[1]]

    if(d==2)
    {
      x[[2]] <- BiCopHinv2(w[[2]],
                           v[[1]][[1]],
                           BiCop(RVM$family[1,1],RVM$par,RVM$par2))#$hinv2
    }
    else
    {
      x[[2]] <- BiCopHinv2(w[[2]],
                        v[[1]][[1]],
                        BiCop(RVM$family[d+1-(1),1],RVM$par[d+1-(1),1],RVM$par2[d+1-(1),1]))#$hinv2
    }
    v[[2]][[1]] <- x[[2]]
  }
  else if(Ncond==1)
  {
    w[[1]] <- cond[[1]]
    for(t in ((Ncond+1):d))
    {
      w[[t]] <- BiCopSim(N,0,0)[,1]
    }

    x[[1]] <-  w[[1]]
    v[[1]][[1]] <- x[[1]]

    if(d==2)
    {
      x[[2]] <- BiCopHinv2(w[[2]],
                          v[[1]][[1]],
                          BiCop(RVM$family[1,1],RVM$par,RVM$par2))#$hinv2
    }
    else
    {
      x[[2]] <- BiCopHinv2(w[[2]],
                          v[[1]][[1]],
                          BiCop(RVM$family[d+1-(1),1],RVM$par[d+1-(1),1],RVM$par2[d+1-(1),1]))#$hinv2
      v[[2]][[1]] <- x[[2]]
    }
  }
  else if(Ncond>=2)
  {
    for(t in 1:Ncond)
    {
      w[[t]] <- cond[[t]]
    }
    for(t in ((Ncond+1):d))
    {
      w[[t]] <- BiCopSim(N,0,0)[,1]
    }
    x[[1]] <-  w[[1]]
    v[[1]][[1]] <- x[[1]]

    x[[2]] <- w[[2]]
    v[[2]][[1]] <- x[[2]]
  }

  if(d==2)
  {
    return(cbind(x[[2]],cond[[1]]))
  }
  else if(d>2)
  {
    v[[2]][[2]] <- BiCopHfunc2(
      v[[1]][[1]],
      v[[2]][[1]],
      RVMFamilySim[d+1-(1),1],###
      RVM$par[d+1-(1),1],
      RVM$par2[d+1-(1),1])#$hfunc2

    for(i in 3:d)
    {
      v[[i]] <- list()
      v[[i]][[1]] <- w[[i]]

      if(i>Ncond)
      {
        for(k in (i-1):2)
        {
          v[[i]][[1]] <- BiCopHinv2(v[[i]][[1]],
                                   v[[i-1]][[2*k-2]],
                                   BiCop(RVM$family[d+1-(k),(i-k)],RVM$par[d+1-(k),(i-k)],RVM$par2[d+1-(k),(i-k)]))#$hinv2
        }

        v[[i]][[1]] <- BiCopHinv2(v[[i]][[1]],
                                 v[[i-1]][[1]],
                                 BiCop(RVM$family[d+1-(1),(i-1)],RVM$par[d+1-(1),(i-1)],RVM$par2[d+1-(1),(i-1)]))#$hinv2
      }

      x[[i]] <- v[[i]][[1]]
      if(i == d)
      {break}
      else
      {
        v[[i]][[2]] <- BiCopHfunc2(
          v[[i-1]][[1]],
          v[[i]][[1]],
          RVMFamilySim[d+1-(1),i-1],###
          RVM$par[d+1-(1),i-1],
          RVM$par2[d+1-(1),i-1])#$hfunc2


        v[[i]][[3]] <- BiCopHfunc2(
          v[[i]][[1]],
          v[[i-1]][[1]],
          RVMFamilySim[d+1-(1),i-1],###
          RVM$par[d+1-(1),i-1],
          RVM$par2[d+1-(1),i-1])#$hfunc2

        if(i>3)
        {
          for(j in 2:(i-2))
          {
            v[[i]][[2*j]] <- BiCopHfunc2(
              v[[i-1]][[2*j-2]],
              v[[i]][[2*j-1]],
              RVMFamilySim[d+1-(j),i-j],###
              RVM$par[d+1-(j),i-j],
              RVM$par2[d+1-(j),i-j])#$hfunc2
            v[[i]][[2*j+1]] <- BiCopHfunc2(
              v[[i]][[2*j-1]],
              v[[i-1]][[2*j-2]],
              RVMFamilySim[d+1-(j),i-j],###
              RVM$par[d+1-(j),i-j],
              RVM$par2[d+1-(j),i-j])#$hfunc2
          }

        }
        v[[i]][[2*i-2]] <- BiCopHfunc2(
          v[[i-1]][[2*i-4]],
          v[[i]][[2*i-3]],
          RVMFamilySim[d+1-(i-1),1],###
          RVM$par[d+1-(i-1),1],
          RVM$par2[d+1-(i-1),1])#$hfunc2
      }
    }

    xORDERED <- x
    matrix0 <- rep(0,d)
    for(k in 1:d)
    {
      matrix0[k] <- RVM$Matrix[k,k]
    }
    for(k in 1:d)
    {
      xORDERED[[k]] <-x[[which(matrix0==k)]]
    }
    x <- xORDERED

    xxx <- matrix(0,nrow=N,ncol=d)
    for(k in 1:d)
    {
      xxx[,k] <- x[[k]]
    }
    return(xxx)
  }
}















PossibleCDVine4Condition <- function(Nx,Ny,type){
  #this routine gives back the first level of all the possible vines that mathematically allow foor building a conditionated model
  x <- seq(1,(Ny+Nx),1)
  #in x the first variable has to be the x variables, the last the predictors
  Perm <- permn(x)

  if(Nx==0)
  {
    Nx=Ny
    Ny=0
  }
  
  GoodTree <- logical(length = factorial(length(x)))
  for(i in 1:factorial(length(x)))
  {
    j <- 0
    repeat
    {
      j <- j+1
      GoodTree[i] <- TRUE
      if((which(x==Perm[[i]][j]) < x[Ny+1]) | (j==Nx))
        #Perm[[i]][j] is not a predictor | finished)
      {
        if((j==Nx) & (which(x==Perm[[i]][j]) >= x[Ny+1]))
        {
          GoodTree[i] <- GoodTree[i]
        }
        else if(which(x==Perm[[i]][j]) < x[Ny+1])
        {
          GoodTree[i] <- FALSE
        }
        break
      }
    }
  }
  aaa <- which(GoodTree==TRUE)
  MatrixGoodTree <- matrix(0,nrow=factorial(Ny)*factorial(Nx),ncol=length(x))
  for(i in 1:length(aaa))
  {
    MatrixGoodTree[i,] <- Perm[[aaa[i]]]
  }
  
  if(type==1 | type=="CVine")
  {
    bbb=MatrixGoodTree[,ncol(MatrixGoodTree):1]#inversion wih respect to Dvine!
    MatrixGoodTree=bbb
    if(Ny!=1)# Ny!=1 -> no duplications!
    {
      aaa=MatrixGoodTree
      d=Nx+Ny
      for(i in 1:length(aaa[,1]))
      {
        if(aaa[i,2]<aaa[i,1])
        {
          MatrixGoodTree[i,2]=aaa[i,1]
          MatrixGoodTree[i,1]=aaa[i,2]
        }
      }
      MatrixGoodTree=unique(MatrixGoodTree)
      MatrixGoodTree
    }
    return(MatrixGoodTree)
  }
  if(type==2 | type=="DVine")
  {
    if(Ny==0 | Nx==0)# if Ny==0 | Nx==0 there are duplications!
    {
      half=dim(MatrixGoodTree)[1]/2
      MatrixGoodTree=MatrixGoodTree[1:half,]
    }
    return(MatrixGoodTree)
  }
}

























WriteMatrixDvine <- function(FirstLevelDvine){
  ddd <- FirstLevelDvine
  MatrixDvine <- matrix(0,nrow=length(ddd),ncol=length(ddd))


  for(i in 1:length(ddd))
  {
    MatrixDvine[i,i] <- ddd[i]
    if(i!=1)
    {
      for(j in 1:(i-1))
      {
        MatrixDvine[i,j] <- ddd[length(ddd)+j-i+1]
      }
    }
  }
  return(MatrixDvine)
}











FirstDVineCond <- function(data,Nx,treecrit="AIC",selectioncrit="AIC",familyset = NA,indeptest = FALSE, level = 0.05, se = FALSE,
                           rotations = TRUE, method = "mle"){
  Ny <- dim(data)[2]-Nx
  if(!is.numeric(familyset)[1])
  {
    familyset=c(1,2,3,4,5,6,7,8,9,10,13,14,16,17,18,19,20,23,24,26,27,28,29,30,33,34,36,37,38,39,
                40,104,114,124,134,
                204,
                214,
                224,
                234
    )
  }

  Dvine <- list()
  PossibleDvineFirstLevel <- PossibleCDVine4Condition(Nx,Ny,"DVine")
  NumbDVine <- length(PossibleDvineFirstLevel[,1])#factorial(Ny)*factorial(Nx)

  
    for(i in 1:NumbDVine)
    {
      print(paste("Number D-vine analysed: ",i,"/",NumbDVine,sep=""))
      Dvine[[i]] <- WriteMatrixDvine(PossibleDvineFirstLevel[i,])

      aaa <- RVineCopSelect(data,familyset,Dvine[[i]],selectioncrit,indeptest = indeptest, level = level, se = se,
                            rotations = rotations, method = method)
      Dvine[[i]] <- aaa
    }

    TestAIC <- rep(1000000,NumbDVine)
    TestBIC <- rep(1000000,NumbDVine)

    for(i in 1:NumbDVine)
    {

      TestAIC[i] <- RVineAIC(data, Dvine[[i]], par = Dvine[[i]]$par, par2 = Dvine[[i]]$par2)$AIC
      TestBIC[i] <- RVineBIC(data, Dvine[[i]], par = Dvine[[i]]$par, par2 = Dvine[[i]]$par2)$BIC
    }
    TestS <- data.frame(tree=seq(1,NumbDVine,1),AIC=TestAIC,BIC=TestBIC)
    BestAicDVine <- Dvine[[which(TestS$AIC==min(TestS$AIC))[1]]]
    BestBicDVine <- Dvine[[which(TestS$BIC==min(TestS$BIC))[1]]]
    BestDVine <- list()

    if(treecrit=="AIC")
    {
      BestDVine[[1]] <- BestAicDVine$Matrix
      BestDVine[[2]] <- min(TestS$AIC)
      if(length(which(TestS$AIC==max(TestS$AIC)))>1)
      {
        print(paste("Attention: there are ", length(which(TestS$AIC==max(TestS$AIC)))," vines (over a total of ", NumbDVine,") equally ranked as the best. Consider to use CDVineCondRank to visulaize the full ranking",sep=""))
      }
    }
    if(treecrit=="BIC")
    {
      BestDVine[[1]] <- BestBicDVine$Matrix
      BestDVine[[2]] <- min(TestS$BIC)
      if(length(which(TestS$BIC==max(TestS$BIC)))>1)
      {
        print(paste("Attention: there are ", length(which(TestS$BIC==max(TestS$BIC)))," vines (over a total of ", NumbDVine,") equally ranked as the best. Consider to use CDVineCondRank to visulaize the full ranking",sep=""))
      }
    }
  return(BestDVine)
}























WriteMatrixCvine <- function(FirstLevelCvine){
  d <- FirstLevelCvine
  MatrixCvine <- matrix(0,nrow=length(d),ncol=length(d))


  for(i in 1:length(d))
  {
    MatrixCvine[i,(1:i)] <- d[i]
  }
  return(MatrixCvine)
}















FirstCVineCond <- function(data,Nx,treecrit="AIC",selectioncrit="AIC",familyset = NA,indeptest = FALSE, level = 0.05, se = FALSE,
                           rotations = TRUE, method = "mle"){
  Ny <- dim(data)[2]-Nx
  if(!is.numeric(familyset)[1])
  {
    familyset=c(1,2,3,4,5,6,7,8,9,10,13,14,16,17,18,19,20,23,24,26,27,28,29,30,33,34,36,37,38,39,
                40,104,114,124,134,
                204,
                214,
                224,
                234
    )
  }
  
  Cvine <- list()
  PossibleCvineFirstLevel <- PossibleCDVine4Condition(Nx,Ny,"CVine")
  NumbCVine <- length(PossibleCvineFirstLevel[,1])


  
    for(i in 1:NumbCVine)
    {
      print(paste("Number C-vine analysed: ",i,"/",NumbCVine,sep=""))
      Cvine[[i]] <- WriteMatrixCvine(PossibleCvineFirstLevel[i,])

      aaa <- RVineCopSelect(data,familyset,Cvine[[i]],selectioncrit,indeptest = indeptest, level = level, se = se,
                            rotations = rotations, method = method)
      Cvine[[i]] <- aaa
    }

    TestAIC <- rep(1000000,NumbCVine)
    TestBIC <- rep(1000000,NumbCVine)

    for(i in 1:NumbCVine)
    {
      TestAIC[i] <- RVineAIC(data, Cvine[[i]], par = Cvine[[i]]$par, par2 = Cvine[[i]]$par2)$AIC
      TestBIC[i] <- RVineBIC(data, Cvine[[i]], par = Cvine[[i]]$par, par2 = Cvine[[i]]$par2)$BIC
    }
    TestS <- data.frame(tree=seq(1,NumbCVine,1),AIC=TestAIC,BIC=TestBIC)
    BestAicCVine <- Cvine[[which(TestS$AIC==min(TestS$AIC))[1]]]
    BestBicCVine <- Cvine[[which(TestS$BIC==min(TestS$BIC))[1]]]
    BestCVine <- list()

    if(treecrit=="AIC")
    {
      BestCVine[[1]] <- BestAicCVine$Matrix
      BestCVine[[2]] <- min(TestS$AIC)
      if(length(which(TestS$AIC==max(TestS$AIC)))>1)
      {
        print(paste("Attention: there are ", length(which(TestS$AIC==max(TestS$AIC)))," vines (over a total of ", NumbCVine,") equally ranked as the best. Consider to use CDVineCondRank to visulaize the full ranking",sep=""))
      }
    }
    if(treecrit=="BIC")
    {
      BestCVine[[1]] <- BestBicCVine$Matrix
      BestCVine[[2]] <- min(TestS$BIC)
      if(length(which(TestS$BIC==max(TestS$BIC)))>1)
      {
        print(paste("Attention: there are ", length(which(TestS$BIC==max(TestS$BIC)))," vines (over a total of ", NumbCVine,") equally ranked as the best. Consider to use CDVineCondRank to visulaize the full ranking",sep=""))
      }
    }
  return(BestCVine)
}
















































RankDVineCond <- function(data,Nx,treecrit="AIC",selectioncrit="AIC",familyset = NA,indeptest = FALSE, level = 0.05, se = FALSE,
                          rotations = TRUE, method = "mle"){
  Ny <- dim(data)[2]-Nx
  if(!is.numeric(familyset)[1])
  {
    familyset=c(1,2,3,4,5,6,7,8,9,10,13,14,16,17,18,19,20,23,24,26,27,28,29,30,33,34,36,37,38,39,
                40,104,114,124,134,
                204,
                214,
                224,
                234
    )
  }

  Dvine <- list()
  PossibleDvineFirstLevel <- PossibleCDVine4Condition(Nx,Ny,"DVine")
  NumbDVine <- length(PossibleDvineFirstLevel[,1])

  
    for(i in 1:NumbDVine)
    {
      print(paste("Number D-vine analysed: ",i,"/",NumbDVine,sep=""))
      Dvine[[i]] <- WriteMatrixDvine(PossibleDvineFirstLevel[i,])

      aaa <- RVineCopSelect(data,familyset,Dvine[[i]],selectioncrit,indeptest = indeptest, level = level, se = se,
                            rotations = rotations, method = method)
      Dvine[[i]] <- aaa
    }

    TestAIC <- rep(1000000,NumbDVine)
    TestBIC <- rep(1000000,NumbDVine)

    for(i in 1:NumbDVine)
    {

      TestAIC[i] <- RVineAIC(data, Dvine[[i]], par = Dvine[[i]]$par, par2 = Dvine[[i]]$par2)$AIC
      TestBIC[i] <- RVineBIC(data, Dvine[[i]], par = Dvine[[i]]$par, par2 = Dvine[[i]]$par2)$BIC
    }

    if(treecrit=="AIC")
    {
      TestS <- data.frame(tree=seq(1,NumbDVine,1),AIC=TestAIC)
      TestS <- TestS[order(TestS$AIC, decreasing =FALSE),]
    }
    if(treecrit=="BIC")
    {
      TestS <- data.frame(tree=seq(1,NumbDVine,1),BIC=TestBIC)
      TestS <- TestS[order(TestS$BIC, decreasing =FALSE),]
    }
    Ranking <- list()

    Ranking[[1]] <- TestS
    Ranking[[2]] <- list()
    for(i in 1:NumbDVine)
    {
      Ranking[[2]][[i]] = Dvine[[TestS$tree[i]]]
    }


    Ranking[[1]]$tree=sort(Ranking[[1]]$tree)
    if(treecrit=="AIC")
    {
      Ranking[[1]]=data.frame(tree=Ranking[[1]]$tree,AIC=Ranking[[1]]$AIC)
    }
    if(treecrit=="BIC")
    {
      Ranking[[1]]=data.frame(tree=Ranking[[1]]$tree,BIC=Ranking[[1]]$BIC)
    }
  return(Ranking)
}








RankCVineCond <- function(data,Nx,treecrit="AIC",selectioncrit="AIC",familyset = NA,indeptest = FALSE, level = 0.05, se = FALSE,
                          rotations = TRUE, method = "mle"){
  Ny <- dim(data)[2]-Nx
  if(!is.numeric(familyset)[1])
  {
    familyset=c(1,2,3,4,5,6,7,8,9,10,13,14,16,17,18,19,20,23,24,26,27,28,29,30,33,34,36,37,38,39,
                40,104,114,124,134,
                204,
                214,
                224,
                234
    )
  }

  Cvine <- list()
  PossibleCvineFirstLevel <- PossibleCDVine4Condition(Nx,Ny,"CVine")
  NumbCVine <- length(PossibleCvineFirstLevel[,1])

  
    for(i in 1:NumbCVine)
    {
      print(paste("Number C-vine analysed: ",i,"/",NumbCVine,sep=""))
      Cvine[[i]] <- WriteMatrixCvine(PossibleCvineFirstLevel[i,])

      aaa <- RVineCopSelect(data,familyset,Cvine[[i]],selectioncrit,indeptest = indeptest, level = level, se = se,
                            rotations = rotations, method = method)
      Cvine[[i]] <- aaa
    }

    TestAIC <- rep(1000000,NumbCVine)
    TestBIC <- rep(1000000,NumbCVine)

    for(i in 1:NumbCVine)
    {

      TestAIC[i] <- RVineAIC(data, Cvine[[i]], par = Cvine[[i]]$par, par2 = Cvine[[i]]$par2)$AIC
      TestBIC[i] <- RVineBIC(data, Cvine[[i]], par = Cvine[[i]]$par, par2 = Cvine[[i]]$par2)$BIC
    }

    if(treecrit=="AIC")
    {
      TestS <- data.frame(tree=seq(1,NumbCVine,1),AIC=TestAIC)
      TestS <- TestS[order(TestS$AIC, decreasing =FALSE),]
    }
    if(treecrit=="BIC")
    {
      TestS <- data.frame(tree=seq(1,NumbCVine,1),BIC=TestBIC)
      TestS <- TestS[order(TestS$BIC, decreasing =FALSE),]
    }
    Ranking <- list()

    Ranking[[1]] <- TestS
    Ranking[[2]] <- list()
    for(i in 1:NumbCVine)
    {
      Ranking[[2]][[i]] = Cvine[[TestS$tree[i]]]
    }

    Ranking[[1]]$tree=sort(Ranking[[1]]$tree)
    if(treecrit=="AIC")
    {
      Ranking[[1]]=data.frame(tree=Ranking[[1]]$tree,AIC=Ranking[[1]]$AIC)
    }
    if(treecrit=="BIC")
    {
      Ranking[[1]]=data.frame(tree=Ranking[[1]]$tree,BIC=Ranking[[1]]$BIC)
    }
  return(Ranking)
}









PossibleCVineMatrixCond <- function(data,Nx){
  Ny <- dim(data)[2]-Nx
  
  Cvine <- list()
  PossibleCvineFirstLevel <- PossibleCDVine4Condition(Nx,Ny,"CVine")
  NumbCVine <- length(PossibleCvineFirstLevel[,1])
  
  for(i in 1:NumbCVine)
  {
    Cvine[[i]] <- WriteMatrixCvine(PossibleCvineFirstLevel[i,])
  }
  return(Cvine)
}




PossibleDVineMatrixCond <- function(data,Nx){
  Ny <- dim(data)[2]-Nx
  
  Dvine <- list()
  PossibleDvineFirstLevel <- PossibleCDVine4Condition(Nx,Ny,"DVine")
  NumbDVine <- length(PossibleDvineFirstLevel[,1])

  for(i in 1:NumbDVine)
  {
    Dvine[[i]] <- WriteMatrixDvine(PossibleDvineFirstLevel[i,])
  }
  return(Dvine)
}








