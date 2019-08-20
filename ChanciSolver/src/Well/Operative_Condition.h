#ifndef OPERATIVE_CONDITION_H
#define OPERATIVE_CONDITION_H

#include <string>

class Operative_Condition{
 private:
    
    std::string _type;
    double _value;
    double _next_change;
 public:
    
    void type(std::string type){_type=type;};
    void value(double value){_value=value;};
    void next_change(double _next_change){_next_change=_next_change;};

    const std::string& type() const {return _type;};
    const double& value() const {return _value;};
    const double& nextChange() const {return _next_change;};
    
};

#endif /* OPERATIVE_CONDITION_H */
