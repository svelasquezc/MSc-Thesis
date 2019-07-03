#ifndef VALUE_READER_H
#define VALUE_READER_H

#include <iostream>
#include <exception>
#include <string>

class Value_Reader{
 public:
    template <class T> void myRead(std::string Message, T& Input, std::string Error){
        while(true){
            std::cout << Message;
            try{
                std::cin >> Input;
                break;
            }catch(std::exception e){
                std::cout << Error <<std::endl;
                continue;
            };
        };
    };
};

#endif /* VALUE_READER_H */
