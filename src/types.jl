"""
VMeasure is type containing the group angle (in radians), group velocity (in m/s, or km/s), error on group angle (in radians), error on log(v).
"""
struct VMeasure
    angle::Real    # group angle (radians)
    value::Real  # group velocity
    err_angle::Real   # error on group angle
    err_value::Real  # error on log(v)
end

"""
    VMeasure(angle::Vector, value::Vector, err_angle::Vector, err_value::Vector)

Return an array of VMeasure, element-wise.
"""
function VMeasure(angle::Vector, value::Vector, err_angle::Vector, err_value::Vector)
    return [VMeasure(angle[i], value[i], err_angle[i], err_value[i]) for i in eachindex(angle)]
end

"""
    VMeasure(angle::Vector, value::Vector, err_angle::Real, err_value::Real)

Return an array of VMeasure, element-wise, using the same errors on angle and log(v).
"""
function VMeasure(angle::Vector, value::Vector, err_angle::Real, err_value::Real)
    return [VMeasure(angle[i], value[i], err_angle, err_value) for i in eachindex(angle)]
end
