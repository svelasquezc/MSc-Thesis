#ifndef VALUE_READER_H
#define VALUE_READER_H

#include <iostream>
#include <exception>
#include <string>
#include <limits>

class Value_Reader{
 public:
    template <class T> static void myRead(std::string Message, T& Input, std::string Error){
        
        std::cout << Message << std::endl;
        while(!(std::cin >> Input)){
        
            std::cin.clear();
            // Throw away the rest of the line
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            std::cout << Error <<std::endl;
        };
                
    };
};

#endif /* VALUE_READER_H */
