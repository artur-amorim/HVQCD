#include <iostream>

using namespace std;

iHQCD <- function(A0 = 5, zmax = 10, h = 0.002, lambda0 = 0.0337462) {
  z     <- seq(exp(-A0), zmax, by = h)
  zmin  <- min(z)
  b0    <- 4.2
  b1    <- 51 * b0^2 /121
  alpha <- 2
  aa    <- (3/8) * (alpha - 1) / alpha
  LAdS  <- 1

  # Solves the Gursoy, Kiritsis and Nitti model.
  # Finds A and lambda as a function of z.
  solveIHQCD  <- function(A0 = 5, zmax = 10, h = 0.002) {
    # save the possible changes
    z     <<- seq(exp(-A0), zmax, by = h)
    zmin  <<- min(z)

    flog.debug(paste('[IHQCD] Solving IHQCD for A0', A0,' in z = [', round(zmin, digits = 4), ',', zmax,';', h, '] grid size', length(z),'\n'))
    pars    <- c(b0 = b0, b1 = b1, aa = aa, LAdS = LAdS)
    initCon <- c(A = A0, lambda = lambda0)

    # differential equations
    BG1 <- function(Time, State, Pars) {
      with(as.list(c(State, Pars)), {
         X <- -b0 * lambda / (3 + 2 * b0 * lambda) - (2 * b0^2 + 3 * b1) * lambda^2 / (9 * (1 + lambda^2) * (1 + (2 * b0^2 + 3 * b1) * log(1 + lambda^2) / (18 * aa)))
         W <- (9/4) * (1 + (2/3) * b0 * lambda)^(2/3) * (1 + (2 * b0^2 + 3 * b1) * log(1 + lambda^2) / (18 * aa))^(4 * aa / 3)
         dA <- -(4/9) * W * exp(A) / LAdS
         dlambda <- -(4/3) * X * W * lambda * exp(A) / LAdS

         return(list(c(dA, dlambda)))
       })
    }

    # solve the problem and then store the relevant variables
    sol          <- ode(func = BG1, y = initCon, parms = pars, times = z)
    A            <- sol[,2]
    lambda       <- sol[,3]

    # compute some extras
    Aspline      <- splinefun(z, A)
    lambdaspline <- splinefun(z, lambda)
    Phi          <- log(lambda)
    Phispline    <- splinefun(z, Phi)
    /* derivatives are useful for different applications, so let's export them as well*/
    Ader1        <- Aspline(z, deriv = 1)
    Ader2        <- Aspline(z, deriv = 2)
    Ader3        <- Aspline(z, deriv = 3)
    lambdader1   <- lambdaspline(z, deriv = 1)
    lambdader2   <- lambdaspline(z, deriv = 2)
    Phider1      <- Phispline(z, deriv = 1)
    Phider2      <- Phispline(z, deriv = 2)
    /* also sometimes the warp factor in the string frame is needed */
    As           <- A + (2/3) * Phi
    AsSpline     <- splinefun(z, As)
    Asder1       <- AsSpline(z, deriv = 1)
    Asder2       <- AsSpline(z, deriv = 2)
    Asder3       <- AsSpline(z, deriv = 3)
    dress        <- exp(2 * A)
    /* compute the tensor glueball potential */
    u2           <- (3/2) * Ader2 + (9/4) * Ader1^2
    LambdaQCD    <- exp(A0 - 1 / (b0 * lambda0)) / (b0 * lambda0)^(b1 / b0^2)
    /* some combinations that appear in the potential */
    aF           <- Phider2
    bF           <- Asder2 - Asder1^2
    cF           <- Phider1^2
    l1_2         <- sqrt(lambda)
    e2As         <- exp(2 * As)
    e2A          <- exp(2 * A)

    /* to avoid z=0 boundary effects the first values are removed
     the biggest problem comes from the derivative of A*/
    len <- length(z)
    st  <- 6
    list(z = z[st:len],
        u2 = u2[st:len],
        A = A[st:len],
        Ader1 = Ader1[st:len],
        Ader2 = Ader2[st:len],
        Ader3 = Ader2[st:len],
        Afun  = Aspline,
        Asder1 = Asder1[st:len],
        Asder2 = Asder2[st:len],
        Asder3 = Asder3[st:len],
        lambda = lambda[st:len],
        lambdader1 = lambdader1[st:len],
        lambdader2 = lambdader2[st:len],
        Phi = Phi[st:len],
        Phider1 = Phider1[st:len],
        Phider2 = Phider2[st:len],
        Phifun  = Phispline,
        As = As[st:len],
        AsSpline = AsSpline,
        dress = dress[st:len],
        aF = aF[st:len],
        bF = bF[st:len],
        cF = cF[st:len],
        l1_2 = l1_2[st:len],
        e2As = e2As[st:len],
        e2A = e2A[st:len],
        h = h,
        A0 = A0,
        zmax = zmax,
        LambdaQCD = LambdaQCD)
  }

}


class IHQCD{
    // Initial value of the warp factor A
    float A0;
    // Initial value of lambda
    float lambda0;
    // The geometry is defined until zmax
    const float zmax;
    // Spacing between adjacent points
    const float h;
    float b0, b1, alpha, aa, LAdS;    
    b0 = 4.2;
    b1 = 51 * b0^2 /121;
    alpha = 2;
    aa  = (3/8) * (alpha - 1) / alpha;
    LAdS = 1;

    public :
        solveIHQCD()
    A0 = 5, zmax = 10, h = 0.002, lambda0 = 0.0337462
};

int main(){
    
    return 0;
}