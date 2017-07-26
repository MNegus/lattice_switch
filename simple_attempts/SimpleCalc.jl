# SimpleCalc.jl
# Alternative to using auto-correlate. Instead reads data at a specified interval

module SimpleCalc

# Function to read in the data file and create an array, removing a certain portion
function trim_data(filename, remove_portion)
  data_array = open(readlines, filename)
  string_arr = data_array[trunc(Int, length(data_array) * remove_portion):end]
  return [parse(Float64, p) for p in string_arr]
end


function calculate(filename, remove_portion; interval_length=100)
  data_arr = trim_data(filename, remove_portion)[1:interval_length:end]
  println(length(data_arr))
  return [mean(data_arr), sqrt(var(data_arr) / length(data_arr))]
end

end