#ifndef PROPERTIES_H
#define PROPERTIES_H
#include <string>
/** \file Properties.h
 * \brief The properties classes header file
 * 
 * This file has a set of classes to wrap primitive data times with get and set methods
 * so that the C# idiom of "properties" can be imititated: note though that calling 
 * the set method will need an extra set of brackets compared to C#
 */
using namespace std;

class IntProperty{
private:
    int value;
public:
    int operator()(){return value;}
    void operator=(int input){value=input;}
};
class StringProperty{
private:
    string value;
public:
    string operator()(){return value;}
    void operator=(string input){value=input;}
};
class UnsignedProperty{
private:
    unsigned value;
public:
    unsigned operator()(){return value;}
    void operator=(unsigned input){value=input;}
};
class FloatProperty{
private:
    float value;
public:
    float operator()(){return value;}
    void operator=(float input){value=input;}
};
class BoolProperty{
private:
    bool value;
public:
    float operator()(){return value;}
    void operator=(bool input){value=input;}
};
class DoubleProperty{
private:
    double value;
public:
    float operator()(){return value;}
    void operator=(double input){value=input;}
};
#endif