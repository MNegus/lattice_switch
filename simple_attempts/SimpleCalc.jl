# SimpleCalc.jl
# Alternative to using auto-correlate. Instead reads data at a specified interval

module SimpleCalc

# Function to read in the data file and create an array, removing a certain portion
function trim_data(data_arr, remove_portion)
  return data_array[trunc(Int, length(data_array) * remove_portion):end]
end


function calculate(data_arr, remove_portion; interval_length=100)
  data_arr = trim_data(data_arr, remove_portion)[1:interval_length:end]
  return [mean(data_arr), sqrt(var(data_arr) / length(data_arr))]
end

end