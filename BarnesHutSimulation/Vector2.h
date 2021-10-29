#pragma once

#define POS_TYPE double

struct Vector2
{
    Vector2() { data[0] = 0; data[1] = 0; }
    Vector2(POS_TYPE x, POS_TYPE y) { data[0] = x; data[1] = y; }
    Vector2(const Vector2& v) { data[0] = v[0]; data[1] = v[1]; }

    POS_TYPE data[2];

    POS_TYPE operator[](int i) const { return data[i]; }
    POS_TYPE& operator[](int i) { return data[i]; }
    Vector2& operator+=(const Vector2& rhs) { data[0] += rhs[0]; data[1] += rhs[1]; return *this; }
    Vector2& operator-=(const Vector2& rhs) { data[0] -= rhs[0]; data[1] -= rhs[1]; return *this; }
    Vector2& operator*=(const float rhs) { data[0] *= rhs; data[1] *= rhs; return *this; }
    Vector2& operator/=(const float rhs) { data[0] /= rhs; data[1] /= rhs; return *this; }
    Vector2 operator+(const Vector2& rhs) { Vector2 v;  v[0] = data[0] + rhs[0]; v[1] = data[1] + rhs[1]; return v; }
    Vector2 operator-(const Vector2& rhs) { Vector2 v;  v[0] = data[0] - rhs[0]; v[1] = data[1] - rhs[1]; return v; }
    Vector2 operator*(const float rhs) { Vector2 v;  v[0] = data[0] * rhs; v[1] = data[1] * rhs; return v; }
    Vector2 operator/(const float rhs) { Vector2 v;  v[0] = data[0] / rhs; v[1] = data[1] / rhs; return v; }

    POS_TYPE lengthSquared() { return data[0] * data[0] + data[1] * data[1]; }
};
