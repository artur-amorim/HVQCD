#include <iostream>
#include <vector>

using namespace std;

template <typename T>
class Stack {
    private:
        vector<T> elems; //elements
    public:
        void push(const T elem);
        void pop();
        T top() const;

        bool empty() const
        {
            return elems.empty() ;
        }

        T dostuff( T x)
        {
            return x ;
        }
};

template <class T>
void Stack<T>::push (const T elem)
{
    elems.push_back(elem);
}

template <class T>
void Stack<T>::pop ()
{
    if (elems.empty())
    {
        throw out_of_range("Stack<>::pop(): empty stack");
    }
    // remove last element
    elems.pop_back();
}

template <class T>
T Stack<T>::top () const
{
    if (elems.empty())
    {
        throw out_of_range("Stack<>::pop(): empty stack");
    }
    return elems.back();
}

int main()
{
    vector<double> vec(5,1);
    vector<double> * p = &vec;
    for(int i = 0; i < 5; i++)
    {
        cout << vec[i] << endl;
    }
    return 1;
}
