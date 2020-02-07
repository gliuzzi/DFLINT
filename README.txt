#######################################
####                               ####
####            DFL int            ####
####                               ####
#### A matlab software to solve    ####
#### black-box problems with all   ####
#### integer variables             ####
####                               ####
#######################################

The provided package solves problems of the form:

    min  f(x)
    s.t. g(x) <= 0

where x is a vector of n unknowns restricted to take integer values,
f and g are black-box functions.

To use the software follow the instruction reported below. 

1) unzip and untar the provided DFLint.tar.gz file in a folder of your choice.
   In the following it is assumed that the folder where DFLint.tar.gz has been 
   unpacked is:  ~/DFLint

2) Within folder ~/DFLint let matlab start

3) define your own matlab function to compute objective and (if any) constraint
   function values. The function MUST return:
   - a single value (the objective function) if no constraint is present,
     see e.g. the file kowalik.m for an example

   - a single value plus a column vector of constraint function values
     (keep in mind that only inequality constraints in the form g(x) <= 0 are
     supported), see e.g. the file davidon2_b.m for an example

4) modify either file example_box.m (if no constraint is present) or
   example_con.m (if constraints are present)

5) at matlab prompt execute

   > example_box

   or

   > example_con

6) If you need further assistance, please do not hesitate and contact Giampaolo Liuzzi
   by sending an email to: giampaolo.liuzzi@iasi.cnr.it
