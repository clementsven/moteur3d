#ifndef STRUCT_H
#define STRUCT_H

struct Vec2i
{
    int x;
    int y;

    Vec2i()
    {
        this->x = 0;
        this->y = 0;
    }

    Vec2i(int x, int y)
    {
        this->x = x;
        this->y = y;
    }

};

struct Vec3f
{
    float x;
    float y;
    float z;

    Vec3f()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }

    Vec3f(float x, float y, float z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

};

#endif // STRUCT_H