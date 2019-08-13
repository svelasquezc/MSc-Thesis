#ifndef EQUATION_H
#define EQUATION_H

#include <memory>
#include <typeinfo>
#include <string>

class Equation_Base{
 protected:    
    int _index;
    bool _status;
    
 public:
    virtual ~Equation_Base() = default; 
    virtual const std::string type() = 0;
    virtual const bool& status() const = 0;
    virtual const int& index() const = 0;
};

template<typename TypeRef>
class Equation : public Equation_Base{
    
 private:    
    std::shared_ptr<TypeRef> _reference;
    
 public:

    typedef TypeRef ReferenceType;

    Equation(std::shared_ptr<TypeRef> reference){
        _reference = reference;
    };
    
    TypeRef& reference() {return *_reference;};

    const std::string type() override { return typeid(ReferenceType).name();};

    const bool& status() const override {return _status;};
    const int& index() const override {return _index;};
};

#endif /* EQUATION_H */
