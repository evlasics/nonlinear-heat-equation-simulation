#pragma once

struct Vec2 {
    double x, y;

    Vec2(double x_ = 0, double y_ = 0)
        : x(x_), y(y_) {}

    Vec2 operator+(const Vec2& o) const {
        return Vec2(x + o.x, y + o.y);
    }

    Vec2 operator-(const Vec2& o) const {
        return Vec2(x - o.x, y - o.y);
    }

    Vec2 operator*(double s) const {
        return Vec2(x * s, y * s);
    }

    Vec2 operator/(double s) const {
        return Vec2(x / s, y / s);
    }
};
