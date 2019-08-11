#ifndef EQUATION_H
#define EQUATION_H

#include <memory>

class Equation_Base{
 protected:    
    int _index;
    bool _status;
    
 public:
    
    virtual const bool& status() const = 0;
    virtual const int& index() const = 0;
};

template<typename TypeRef>
class Equation : Equation_Base{
    
 private:    
    std::shared_ptr<TypeRef> _reference;
    
 public:

    typedef TypeRef ReferenceType;

    Equation(std::shared_ptr<TypeRef> reference){
        _reference = reference;
    };
    
    TypeRef& reference() {return *_reference;};

    virtual const bool& status() const{return _status;};
    virtual const int& index() const{return _index;};
};

#endif /* EQUATION_H */
