# test 
# https://stackoverflow.com/questions/57856246/t-test-with-bootstrap-in-r
# Here's an example of using that function with simulated data 
# where you'd expect a p-value close to 1.

set.seed(0)
df <- data.frame(gender = sample(c('M', 'F'), size=50, replace=T),
                 measure = runif(n=50))

boot.t.test(df[df$gender=='M', 'measure'], df[df$gender=='F', 'measure'], reps=1000)
