# extra functions for gammal the analysis GLMs #######

devresids <- function( m )
{
  # pointwise deviances (squared deviance residuals)
  y <- m$y
  mu <- m$family$linkinv( m$linear.predictors )
  wt <- m$weights
  -2 * wt * ( log( ifelse( y == 0, 1, y / mu ) ) - ( y - mu ) / mu )
  # sign( y - mu ) * sqrt( devresids( m ) ) - residuals( object = m, type = 'deviance' )
}
devres <- function( m )
{
  # -2 * sum( log( y / mu ) - ( y - mu ) / mu )
  sum( devresids( m ) )
}
glogLik <- function( m, dispersion )
{
  # log likelihood
  if ( !( 'glm' %in% class( m ) ) || ( m$family$family != 'Gamma' ) )
  {
    stop( 'please input a gamma glm' )
  }
  
  if ( missing( dispersion ) )
  {
    # use MLE dispersion estimate provided in package MASS by default
    require( MASS )
    phi <- 1 / gamma.shape( m )$alpha
  } else
  {
    # use dispersion estimate provided by the user
    phi <- dispersion
  }
  
  y <- m$y
  mu <- m$family$linkinv( m$linear.predictors )
  a <- 1 / phi
  b <- mu / a
  
  sum( dgamma( x = y, shape = a, scale = b, log = TRUE ) )
  # sum( - a * log( b ) - lgamma( a ) + ( a - 1 ) * log( y ) - y / b  )
  # sum( glogLikp( m, dispersion )
}
glogLikp <- function( m, dispersion )
{
  # pointwise log likelihood
  if ( !( 'glm' %in% class( m ) ) || ( m$family$family != 'Gamma' ) )
  {
    stop( 'please input a gamma glm' )
  }
  
  if ( missing( dispersion ) )
  {
    # use MLE dispersion estimate provided in package MASS by default
    require( MASS )
    phi <- 1 / gamma.shape( m )$alpha
  } else
  {
    # use dispersion estimate provided by the user
    phi <- dispersion
  }
  
  y <- m$y
  mu <- m$family$linkinv( m$linear.predictors )
  a <- 1 / phi
  b <- mu / a
  
  dgamma( x = y, shape = a, scale = b, log = TRUE )
  # - a * log( b ) - lgamma( a ) + ( a - 1 ) * log( y ) - y / b 
}
gAIC <- function( m, dispersion )
{
  # AIC for a gamma GLM with the option to use a custom dispersion, phi
  if ( !( 'glm' %in% class( m ) ) || ( m$family$family != 'Gamma' ) )
  {
    stop( 'please input a gamma glm' )
  }
  
  if ( missing( dispersion ) )
  {
    # use MLE dispersion estimate provided in package MASS
    require( MASS )
    phi <- 1 / gamma.shape( m )$alpha
  } else
  {
    # use dispersion estimate provided by the user
    phi <- dispersion
  }
  
  k <- m$rank + 1 # number of parameters (the extra one is for the shape)
  2.0 * ( k - glogLik( m, phi ) )
}
gAICc <- function( m, dispersion )
{
  # AICc for a gamma GLM with the option to use a custom dispersion, phi
  if ( !( 'glm' %in% class( m ) ) || ( m$family$family != 'Gamma' ) )
  {
    stop( 'please input a gamma glm' )
  }
  
  if ( missing( dispersion ) )
  {
    # use MLE dispersion estimate provided in package MASS
    require( MASS )
    phi <- 1 / gamma.shape( m )$alpha
  } else
  {
    # use dispersion estimate provided by the user
    phi <- dispersion
  }
  
  n <- length( m$y )
  k <- m$rank + 1 # number of parameters (the extra one is for the shape)
  2.0 * ( k * n / ( n - k - 1.0 ) - glogLik( m, phi ) )
}
AICc <- function( m )
{
  # AICc for a gamma GLM with the option to use a custom dispersion, phi
  if ( !( 'glm' %in% class( m ) ) || ( m$family$family != 'Gamma' ) )
  {
    stop( 'please input a gamma glm' )
  }
  n <- length( m$y )
  k <- m$rank + 1 # number of parameters (the extra one is for the shape)
  2.0 * ( k * n / ( n - k - 1.0 ) - as.numeric( logLik( m ) ) )
}
gNagelkerke <- function( m, dispersion = NULL, dispersion.null = NULL )
{
  if ( !( 'glm' %in% class( m ) ) || ( m$family$family != 'Gamma' ) )
  {
    stop( 'please input a gamma glm' )
  }
  
  if ( missing( dispersion ) )
  {
    # use MLE dispersion estimate provided in package MASS
    require( MASS )
    phi <- 1 / gamma.shape( m )$alpha
  } else
  {
    # use dispersion estimate provided by the user
    phi <- dispersion
  }
  m.null <- glm( formula = m$y ~ 1, family = family( m ) )
  if ( missing( dispersion.null ) )
  {
    # use MLE dispersion estimate provided in package MASS
    require( MASS )
    phi.null <- 1 / gamma.shape( m.null )$alpha
  } else
  {
    # use dispersion estimate provided by the user
    phi.null <- dispersion.null
  }
  n <- length( m$y )
  ( 1.0 - exp( - 2.0 * as.numeric( glogLik( m, phi ) - glogLik( m.null, phi.null ) ) / n ) ) / as.numeric( 1.0 - exp( 2.0 * glogLik( m.null, phi.null ) / n ) )
}

#S-R
fsr <- function( xst, xt, sr )
{
  # residual time series after fitting a SR model
  xout <- rep( x = NA, times = length( xst ) )
  if ( sr == 'Ricker' )
  {
    xout <- resid( lm( log( xst ) ~ xt, na.action = na.exclude ) )
    
  } else if ( sr == 'Beverton-Holt' )
  {
    # starting values
    # m <- lm( I( xst * xt ) ~ 1 + xt, na.action = na.exclude )
    
    xout <- resid( nls( log( xst ) ~ a - log( 1 + b * xt ), 
                        algorithm = 'port',
                        lower = c( 0, 0 ), 
                        # start = c( a = 0.5 *log( coef( m )[ 1 ] ), b = 1 / max( xt, na.rm = T ) ),
                        start = c( a = 1.0, b = 1 / max( xt, na.rm = T ) ),
                        na.action = na.exclude, control = nls.control( maxiter = 100, warnOnly = TRUE ) ) )
  } else
  {
    stop( paste( sr, ' is an invalid value for \'sr\'; only \'Ricker\' or \'Beverton-Holt\' curves are available', sep = '' ) )
  }
  
  return( xout )
}

#R2
fdevexpl <- function( m )
{
  # fraction of deviance explained
  m.null <- glm( formula = m$y ~ 1, family = family( m ), data = m$data )
  1.0 - sum( devresids( m ) ) / sum( devresids( m.null ) ) # fraction of deviance explained
}

fgNagelkerke <- function( m )
{
  # Nagelkerke r2
  m.null <- glm( formula = m$y ~ 1, family = family( m ), data = m$data )
  k <- m$df.null + 1
  # ( 1.0 - exp( - 2.0 * as.numeric( sum( glogLikp( m ) ) - sum( glogLikp( m.null ) ) ) / k ) ) /
  #     as.numeric( 1.0 - exp( 2.0 * sum( glogLikp( m.null ) ) / k ) )
  # -expm1( - 2.0 * as.numeric( sum( glogLikp( m ) ) - sum( glogLikp( m.null ) ) ) / k ) /
  #   -expm1( 2.0 * sum( glogLikp( m.null ) ) / k ) 
  -expm1( - 2.0 * sum( glogLikp( m ) - glogLikp( m.null ) ) / k ) /
    -expm1( 2.0 * sum( glogLikp( m.null ) ) / k )
}


fdevexplSSet <- function( m, idx )
{
  # fraction of deviance explained for a subset of the data
  if ( missing( idx ) )
    idx <- rep( T, length( m$y ) )
  m.null <- glm( formula = m$y ~ 1, family = family( m ), data = m$data )
  1.0 - sum( devresids( m )[ idx ] ) / sum( devresids( m.null )[ idx ] ) # fraction of deviance explained
}

fgNagelkerkeSSet <- function( m, idx )
{
  # Nagelkerke r2 for a subset of the data
  if ( missing( idx ) )
    idx <- rep( T, length( m$y ) )
  m.null <- glm( formula = m$y ~ 1, family = family( m ), data = m$data )
  ( 1.0 - exp( - 2.0 * as.numeric( sum( glogLikp( m )[ idx ] ) - sum( glogLikp( m.null )[ idx ] ) ) / sum( idx ) ) ) /
      as.numeric( 1.0 - exp( 2.0 * sum( glogLikp( m.null )[ idx ] ) / sum( idx ) ) )
}


