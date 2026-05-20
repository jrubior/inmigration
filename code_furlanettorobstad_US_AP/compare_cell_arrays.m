% Compact comparison
function identical = compare_cell_arrays(theta, theta2, tolerance)
    if nargin < 3
        tolerance = 1e-12;
    end
    
    if length(theta) ~= length(theta2)
        identical = false;
        return;
    end
    
    identical = true;
    for i = 1:length(theta)
        if isnumeric(theta{i}) && isnumeric(theta2{i})
            if ~isequal(size(theta{i}), size(theta2{i})) || ...
               max(abs(theta{i}(:) - theta2{i}(:))) > tolerance
                identical = false;
                return;
            end
        else
            if ~isequal(theta{i}, theta2{i})
                identical = false;
                return;
            end
        end
    end
end

