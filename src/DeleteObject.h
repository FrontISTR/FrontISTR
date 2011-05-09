//
// DeleteObject.h
//
namespace pmw{
struct DeleteObject{
    template<typename T>
    void operator()(const T* ptr) const
    {
        delete ptr;
    }
};
}
