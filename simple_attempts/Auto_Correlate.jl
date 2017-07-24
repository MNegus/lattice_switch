# auto_correlate.jl
# Performs auto-correlation on data taken from Monte-Carlo switching
module Auto_Correlate

# Function to read in the data file and create an array, removing a certain portion
function trim_data(filename, remove_portion)
  data_array = open(readlines, filename)
  string_arr = data_array[trunc(Int, length(data_array) / 2):end]
  return [parse(Float64, p) for p in string_arr]
end

# Perform auto-correlation to produce an array of uncorrelated data
function autocorrelate(data_array, R_cutoff; interval_length=1)
  mean_data = mean(data_array)
  N = length(data_array)
  denomimator = sum([(value - mean_data)^2 for value in data_array])

  period = 0

  for k = [1:interval_length:N;]
    # println(k)
    R = sum([(data_array[j] - mean_data) * (data_array[j + k] - mean_data) for j = [1:N - k;]]) / denomimator
    # println(R)
    if abs(R) <= R_cutoff
      period = k
      # println(period)
      break
    end
  end

  return data_array[1:period:end]
end

function calculate(filename, remove_portion, R_cutoff; interval_length=1)
  data_arr = trim_data(filename, remove_portion)
  uncorr_data = autocorrelate(data_arr, R_cutoff, interval_length=interval_length)
  return [mean(uncorr_data), sqrt(var(uncorr_data) / length(uncorr_data))]
end

end
