#ifndef __arrays_h__
#define __arrays_h__

#include <cstdint>


#include "cuda_defines.h"



namespace hw
{
    //simulation of fortran arrays
    
    template < typename t, int32_t width, int32_t height >
    class array_2d_mn
    {
        t* m_data;

    public:
        __host__ __device__ explicit array_2d_mn(t* memory_buffer) : m_data( memory_buffer)
        {

        }

        __host__ __device__ inline t& operator() (int32_t col, int32_t row)
        {
            return *( m_data + ( row - 1) * width + ( col - 1 ) );
        }

        __host__ __device__ inline t operator() (int32_t col, int32_t row) const
        {
            return *(m_data + (row - 1) * width + (col - 1));
        }

        __host__ __device__ inline const t* begin() const
        {
            return m_data;
        }

        __host__ __device__ inline t* begin()
        {
            return m_data;
        }

        __host__ __device__ inline const t* end() const
        {
            return begin() + (height * width);
        }

        __host__ __device__ inline t* end()
        {
            return begin() + (height * width);
        }
    };
    
    /*
    template < typename t, int32_t width, int32_t height >
    class array_2d_mn
    {
        t m_data[height][width];
  
    public:
        explicit array_2d_mn(t* memory_buffer)
        {

        }

        inline t& operator() (int32_t col, int32_t row)
        {
            return m_data[row - 1][col - 1];
        }

        inline t operator() (int32_t col, int32_t row) const
        {
           return m_data[row - 1][col - 1];
        }

        inline const t* begin() const
        {
            return &m_data[0][0];
        }

        inline t* begin()
        {
            return &m_data[0][0];
        }

        inline const t* end() const
        {
            return begin() + sizeof(m_data) / sizeof(t);
        }

        inline t* end()
        {
            return begin() + sizeof(m_data) / sizeof(t);
        }
    };
    */


    template < typename t, int32_t size >
    class array_1d_mn
    {
        t* m_data;

    public:

        array_1d_mn(t* memory_buffer) : m_data(memory_buffer)
        {

        }

        inline t& operator() (int32_t index)
        {
            return m_data + index - 1;
        }

        inline t operator() (int32_t index) const
        {
            return m_data + index - 1;
        }

        inline const t* begin() const
        {
            return m_data;
        }

        inline t* begin()
        {
            return m_data;
        }

        inline const t* end() const
        {
            return begin() + size;
        }

        inline t* end()
        {
            return begin() + size;
        }
    };

    
    template < typename t, int32_t size >
    class array_1d_mn_cpu
    {
        t m_data[size];

    public:

        inline t& operator() (int32_t index)
        {
            return m_data[index - 1];
        }

        inline t operator() (int32_t index) const
        {
            return m_data[index - 1];
        }

        inline const t* begin() const
        {
            return &m_data[0];
        }

        inline t* begin()
        {
            return &m_data[0];
        }

        inline const t* end() const
        {
            return begin() + sizeof(m_data) / sizeof(t);
        }

        inline t* end()
        {
            return begin() + sizeof(m_data) / sizeof(t);
        }
    };

    
}


#endif