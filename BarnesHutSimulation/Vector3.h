#pragma once

#include "Vector2.h"

struct Vector3
{
    Vector3() { data[0] = 0; data[1] = 0; data[2] = 0; }
    Vector3(POS_TYPE x, POS_TYPE y, POS_TYPE z) { data[0] = x; data[1] = y; data[2] = z; }
    Vector3(const Vector3& v) { data[0] = v[0]; data[1] = v[1]; data[2] = v[2]; }

    POS_TYPE data[3];

    POS_TYPE operator[](int i) const { return data[i]; }
    POS_TYPE& operator[](int i) { return data[i]; }
    Vector3& operator+=(const Vector3& rhs) { data[0] += rhs[0]; data[1] += rhs[1]; data[2] += rhs[2]; return *this; }
    Vector3& operator-=(const Vector3& rhs) { data[0] -= rhs[0]; data[1] -= rhs[1]; data[2] -= rhs[2]; return *this; }
    Vector3& operator*=(const float rhs) { data[0] *= rhs; data[1] *= rhs; data[2] *= rhs; return *this; }
    Vector3& operator/=(const float rhs) { data[0] /= rhs; data[1] /= rhs; data[2] /= rhs; return *this; }
    Vector3 operator+(const Vector3& rhs) { Vector3 v;  v[0] = data[0] + rhs[0]; v[1] = data[1] + rhs[1]; v[2] = data[2] + rhs[2]; return v; }
    Vector3 operator-(const Vector3& rhs) { Vector3 v;  v[0] = data[0] - rhs[0]; v[1] = data[1] - rhs[1]; v[2] = data[2] - rhs[2]; return v; }
    Vector3 operator*(const float rhs) { Vector3 v;  v[0] = data[0] * rhs; v[1] = data[1] * rhs; v[2] = data[2] * rhs; return v; }
    Vector3 operator/(const float rhs) { Vector3 v;  v[0] = data[0] / rhs; v[1] = data[1] / rhs; v[2] = data[2] / rhs; return v; }

    POS_TYPE lengthSquared() { return data[0] * data[0] + data[1] * data[1] + data[2] * data[2]; }
    POS_TYPE dot(Vector3 v) { return  data[0] * v[0] + data[1] * v[1] + data[2] * v[2]; }
};