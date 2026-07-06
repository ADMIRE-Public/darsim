function result = compare_signs(a, b, tol)
    % returns 1 if a and b have the same sign
    % returns 0 otherwise  
   result = (a >= -tol  & b >= -tol) | (a <= tol & b <= tol);
  end