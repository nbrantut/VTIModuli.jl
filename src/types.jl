struct VMeasure
    angle::Float64    # group angle (radians)
    value::Float64  # group velocity
    err_angle::Float64   # error on group angle
    err_value::Float64  # error on log(v)
end
