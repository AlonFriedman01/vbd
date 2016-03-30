#' This function allows you to solve a simple vector borne diseases mathematical model
#' with delay (SIRSvIv). To try the solver please run the following example.
#'
#' Parameters and initial conditions.
#' @param
#'  tau =  Time delay
#' @param N = Total host population
#' @param beta = Host infection rate
#' @param r =  Host recovery rate
#' @param lambda = Number of new suceptible vectors per time
#' @param alpha = Vector infection rate
#' @param mu = Vector death rate
#' @examples
#' parameters <- c(mu = 1 / (70 * 365), beta = 520 / 365,
#'                 lambda = 1 / 14, alpha = 1 / 7, tau =1, N = 500, r = 1/6 )
#' initials <- c(S = 500,  I = 0, R = 0, Sv =0 , Iv = 1)
#' # Solve and plot.
#' S<- SIRSvIv(pars = parameters, init = initials, time = seq(from = 0, to = 10, by =.5))
#' plot(S$results)
#' @export

SIRSvIv <- function(pars = NULL, init = NULL, time = NULL, ...) {
  if (is.null(pars)) {
    stop("undefined 'pars'")
  }
  if (is.null(pars)) {
    stop("undefined 'inits'")
  }
  if (is.null(pars)) {
    stop("undefined 'time'")
  }
  function1 <- function(pars = NULL, init = NULL, time = NULL) {
    function2 <- function(time, init, pars) {
      with(as.list(c(init, pars)), {
        if(time>=tau)
          lag1 <- lagvalue(time - tau)
        else
          lag1 <- c(1,1,1,1,1)
        dS  = -beta*Iv * S/N
        dI  = beta*lag1[1]*lag1[5]/N - r*I
        dR  = r*I
        dSv = lambda - alpha*I*Sv/N - mu*Sv
        dIv = alpha*I*Sv/N - mu*Iv
        list(c(dS, dI, dR, dSv, dIv))
      })
    }
    init <- c(init['S'],init['I'], init['R'], init['Sv'],init['Iv'])
    output <- dede(times = time,
                  func = function2,
                  y = init, parms = pars, ...)
    return(output)
  }

  output <- function1(pars = pars, init = init, time = time)

  return(list(model = function1,
              pars = pars,
              init = init,
              time = time,
              results = as.data.frame(output)))

}

