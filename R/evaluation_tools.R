missing_arg <- function() quote(expr=)

evaluate_single_param <- function(samples, param_name, indices, trueValue) {
  if(is.null(indices)) {
    param_samples = samples[[param_name]];
  }
  else {
    #A magic form to run samples[[param_name]][,indices[1], ... , indices[N]] based on the length of indices
    param_samples = do.call(`[`,append(list(samples[[param_name]],missing_arg()),indices))
  }

  indices_str = do.call(paste, append(indices, list(sep = ",")))
  fullName =   paste0(param_name,"[", indices_str, "]")


  return(data.frame(
    param_name = fullName,
    trueValue = trueValue,
    median = median(param_samples),
    IQR = IQR(param_samples),
    quantile = ecdf(param_samples)(trueValue)
  ))
}


evaluate_all_params <- function(samples, true_params) {
  next_element = 1
  result = list();
  for(param_name in names(true_params)) {
    if(!param_name %in% names(samples)) {
      next;
    }
    param_values = get(param_name, true_params);
    if(is.null(dim(param_values)) || (length(dim(param_values)) == 1)) {
      if(length(param_values) > 1) {
        for(i in 1:length(param_values)) {
          result[[next_element]] = evaluate_single_param(samples, param_name, list(i), param_values[i])
          next_element = next_element + 1
        }
      } else {
        result[[next_element]] = evaluate_single_param(samples, param_name, NULL, param_values)
        next_element = next_element + 1
      }
    }
    else {
      if(length(dim(param_values)) == 2)
      {
        for(i in 1:dim(param_values)[1]) {
          for(j in 1:dim(param_values)[2]) {
            result[[next_element]] = evaluate_single_param(samples, param_name, list(i,j), param_values[i,j])
            next_element = next_element + 1
          }
        }
      } else {
        stop("3+ dimensional parameters not supported yet");
      }
    }
  }
  return(do.call(rbind.data.frame, result));
}

evaluation_summary <- function(samples, true_params, printParamsResults = TRUE) {
  eval_result = evaluate_all_params(samples, true_params);
  if(printParamsResults) {
    print(eval_result);
  }
  quantiles = eval_result$quantile;
  within25 = mean(quantiles >= 0.375 & quantiles <= 0.625);
  within50 = mean(quantiles >= 0.25 & quantiles <= 0.75);
  within95 = mean(quantiles >= 0.025 & quantiles <= 0.975);
  cat("\nWithin 25% interval:", within25,"\nWithin 50% interval:", within50, "\nWithin 95% interval:",within95,"\n")
}
