-- Simple statistics :: Jamie Hudspith Fri. 21 Sep 2012

-- Creates an average of the data
function calc_average( data )
   local sum = 0
   local n = 0
   for i in ipairs( data ) do
      sum = sum + data[ i ]
      n = n + 1
   end
   return sum / n
end

-- Calculates the standard error of the data
function calc_standerr( data )
   local xbar = calc_average( data ) 
   local var = 0 
   local n = 0
   for i in ipairs( data ) do
      local val = data[i] - xbar 
      var = var + ( val * val )
      n = n + 1
   end
   return xbar , math.sqrt( var ) / n
end

-- Calculates the variance of the data
function calc_variance( data )
   local xbar = calc_average( data ) 
   local var = 0 
   local n = 0
   for i in ipairs( data ) do
      local val = data[i] - xbar 
      var = var + ( val * val )
      n = n + 1
   end
   return xbar , var / ( n - 1 ) 
end

-- Calculates the median of the data
function calc_median( data )
   sorted_data = data
   table.sort( sorted_data )
   local n = #sorted_data
   if n%2 == 1 then
      return sorted_data[(n+1)/2] , sorted_data[ 1 ]
   else
      return ( sorted_data[ n/2 ] + sorted_data[ n/2 + 1] )/2.0
   end
end

-- Calculates the median of the data
function calc_maxmin( data )
   sorted_data = data
   table.sort( sorted_data )
   local n = #sorted_data
   return sorted_data[1] , sorted_data[n]
end

-- Creates a histogram of the raw data
function histogram ( data , min , max , nob )

   local histo = {} -- what we return
   local binsize = ( max - min ) / nob -- size of each bin
   local binmids = {}

   -- loop the number of bins
   for bin = 0 , nob do
      
      -- set the upper bound --
      upbound = min + bin * binsize
      local midpoint = upbound - binsize / 2 
      binmids[ bin ] = midpoint ;
      histo[ bin ] = 0

      for i in ipairs( data ) do

         local val = tonumber( data[ i ] )

	 if val < upbound then
	    if val > ( upbound - binsize ) then
	       --print( values )
	       histo[bin] = histo[bin] + 1
	    end
	 end
	 
      end	
   end
   
   -- print the histogram
   local median = 0
   for l in ipairs( histo ) do
      print( binmids[l] , histo[l] )
      local val = tonumber( histo[l] )
      local valcomp = tonumber( histo[ median ] )
      if val > valcomp then median = l end
   end
   print("HISTOGRAM MODE :: ", binmids[ median ] , binsize / 2 )

   -- return the histogram
   return binmids , histo ;
end

-- rounding function needed for the index
function round (x)
   if x >= 0 then
      return math.floor (x + 0.5)
   end  -- if positive
   return math.ceil (x - 0.5)
end -- function round

-- estimate error from distribution
function distrib_err( data , boots , confidence )
   -- sort the data
   local sorted = data
   table.sort( sorted )
   local omitted = round( ( 100 - confidence ) / 2.0 )
   local bottom  = round( ( omitted * boots / 100 ) )
   local top     = round( ( boots - 1 - bottom ) )
   return 0.5 * ( sorted[top] - sorted[bottom] ) 
end

-- Jackknife the data
function jackknife_data( data ) 

   -- initial setup
   local jackdata = {}
   local average = 0 
   local n = 0

   -- loop the data
   for i in ipairs( data ) do
      local count = 0
      jackdata[ i ] = 0

      -- loop through the data performing some sort of average
      for j in ipairs( data ) do
	 -- first of all the average
	 if i ~= j then
	    jackdata[ i ] = jackdata[ i ] + data[ j ]
	    count = count + 1
	 end

      end
      jackdata[ i ] = jackdata[ i ] / count

      average = average + jackdata[ i ]
      n = n + 1
   end
   average = average / n  

   -- compute the error
   local err_jackdata = 0
   local factor = 1.0 - 1.0 / n
   for i in ipairs( jackdata ) do	
      local cache = ( jackdata[i] - average )
      err_jackdata = err_jackdata + factor * cache * cache 
   end

   -- return the average and the error like the variance function
   return average , math.sqrt(err_jackdata)
end

-- bootstrap the data ...
function bootstrap( data , confidence , boots , bootdata )
   -- initial setup
   local average = 0 
   local n = 0

   -- get the data size ...
   for j in ipairs( data ) do	
      average = average + data[ j ]
      n = n + 1
   end
   average = average / n  

   -- initialise the seed, could be better
   math.randomseed( 12345 )    

   -- loop the bootstraps	
   for i = 0 , boots do

      -- initialise the bootstrapped data to 0
      bootdata[ i ] = 0

      -- loop through the data performing some sort of average
      for j in ipairs( data ) do

	 -- first of all the average of this bootstrap
	 local idx = math.random( n )
	 --print( idx )
	 bootdata[ i ] = bootdata[ i ] + data[ idx ]
	 
      end

      -- average this boot
      bootdata[ i ] = bootdata[ i ] / n

   end

   -- bootstrapped average and variance
   local av, var = calc_variance( bootdata ) 

   -- return the average and the symmetric error
   return av, distrib_err( bootdata , boots , confidence )
end
