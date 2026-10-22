/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025-2025 OpenCFD Ltd.
    Copyright (C) 2026 Keysight Technologies
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Simple compilation tests and access for UPstream types

\*---------------------------------------------------------------------------*/

#include "pTraits.H"
#include "contiguous.H"
#include "FixedList.H"
#include "boolVector.H"  // A FixedList pretending to be a vector
#include "barycentric.H"
#include "complex.H"
#include "vector.H"
#include "vector2D.H"
#include "tensor.H"
#include "treeBoundBox.H"
#include "uLabel.H"
#include "MinMax.H"
#include "Switch.H"
#include "Tuple2.H"
#include "IOstreams.H"
#include "UPstream.H"

#include <functional>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Add in some extras from functional

//- Map std::plus to \c UPstream::opCodes::op_sum
template<>
struct UPstream_opType<std::plus<void>> : std::true_type
{
    static constexpr auto opcode_id = UPstream::opCodes::op_sum;
};


//- Map 'signed char' to UPstream::dataTypes::type_byte
//  Caution with: may be identical to int8_t mapping!!
#if 0
template<>
struct UPstream_alias_dataType<signed char> : std::true_type
{
    using base = char;
    static constexpr auto datatype_id = UPstream::dataTypes::type_byte;
};
#endif

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Just for debugging
static const Foam::List<std::string> dataType_names
({
    "byte",
    "int16",
    "int32",
    "int64",
    "uint16",
    "uint32",
    "uint64",
    "float",
    "double",
    "long_double",

    "int_int",
    "long_int",
    "float_int",
    "double_int",

    "float[3]",
    "double[3]",
    "float[6]",
    "double[6]",
    "float[9]",
    "double[9]"
});

// Just for debugging
static const Foam::List<std::string> opType_names
({
    "op_min",
    "op_max",
    "op_sum",
    "op_prod",
    "op_bool_and",
    "op_bool_or",
    "op_bool_xor",
    "op_bit_and",
    "op_bit_or",
    "op_bit_xor",
    "op_minloc",
    "op_maxloc",
    "op_replace",
    "op_no_op"
});


using namespace Foam;

void printDataTypeId(UPstream::dataTypes datatype_id)
{
    if (datatype_id != UPstream::dataTypes::invalid)
    {
        const int index = int(datatype_id);
        if (index < dataType_names.size())
        {
            Info<< dataType_names[index];
        }
        else
        {
            Info<< '(' << index << ')';
        }
    }
}


void printOpCodeId(UPstream::opCodes opcode_id)
{
    if (opcode_id != UPstream::opCodes::invalid)
    {
        const int index = int(opcode_id);
        if (index < opType_names.size())
        {
            Info<< ':' << opType_names[index].c_str();
        }
        else
        {
            Info<< '(' << index << ')';
        }
    }
    else
    {
        Info<< "(null)";
    }
}


template<class T, bool showSize = false>
void printTypeName()
{
    // Both float and double have pTraits typeName = "scalar"!
    if constexpr (std::is_same_v<float, std::remove_cv_t<T>>)
    {
        Info<< "<float>";
    }
    else if constexpr (std::is_same_v<double, std::remove_cv_t<T>>)
    {
        Info<< "<double>";
    }
    else if constexpr (pTraits_has_typeName<T>)
    {
        Info<< pTraits<std::remove_cv_t<T>>::typeName;
    }
    else
    {
        Info<< error::demangle<T>();
    }
    if constexpr (showSize)
    {
        Info<< " (" << sizeof(T) << " bytes)";
    }
}


template<class Type, bool UseTypeName = true>
void printPstreamTraits(std::string_view name = std::string_view())
{
    Info<< "========" << nl;
    Info<< "type: ";
    if (!name.empty())
    {
        Info<< name << ' ';
    }

    if constexpr (UseTypeName)
    {
        printTypeName<Type, true>();
    }
    else
    {
        Info<< error::demangle<Type>()
            << " (" << sizeof(Type) << " bytes)";
    }

    {
        using cmpt = typename Foam::pTraits_cmptType<Type>::type;

        if constexpr (!std::is_same_v<Type, cmpt>)
        {
            Info<< ", cmpt:";

            if constexpr (UseTypeName)
            {
                printTypeName<cmpt, true>();
            }
            else
            {
                Info<< error::demangle(typeid(cmpt))
                    << " (" << sizeof(cmpt) << " bytes)";
            }
        }
    }


    Info<< nl
        << "  is_contiguous:"
        << is_contiguous<Type>::value;

    if constexpr (UPstream_mpi_dataType<Type>::value)
    {
        Info<< ", is_mpi=("
            << int(UPstream_mpi_dataType<Type>::datatype_id) << ')';
    }
    else
    {
        std::cout << ", is_mpi=(null)";
    }
    if constexpr (UPstream_user_dataType<Type>::value)
    {
        Info<< ", is_user=("
            << int(UPstream_user_dataType<Type>::datatype_id) << ')';
    }
    else
    {
        std::cout << ", is_user=(null)";
    }
    if constexpr (UPstream_any_dataType<Type>::value)
    {
        Info<< ", is_any=("
            << int(UPstream_any_dataType<Type>::datatype_id) << ')';
    }
    else
    {
        std::cout << ", is_any=(null)";
    }

    // Any aliases?
    if constexpr
    (
        UPstream_alias_dataType<Type>::value
     && !UPstream_mpi_dataType<Type>::value
    )
    {
        Info<< ", alias=("
            << int(UPstream_alias_dataType<Type>::datatype_id) << ')';
    }

    Info<< " base-type:" << int(UPstream_basic_dataType<Type>::datatype_id)
        << " data-type:" << int(UPstream_dataType<Type>::datatype_id)
        << nl;

    if constexpr (UPstream_basic_dataType<Type>::value)
    {
        Info<< " base-type=";
        printDataTypeId(UPstream_basic_dataType<Type>::datatype_id);
    }
    else if constexpr (UPstream_dataType<Type>::value)
    {
        Info<< " data-type=";
        printDataTypeId(UPstream_dataType<Type>::datatype_id);
    }

    {
        // Use element or component type (or byte-wise) for data type
        using base = typename UPstream_dataType<Type>::base;

        // The sizing factor is constexpr
        constexpr std::streamsize count = UPstream_dataType<Type>::size(1);

        Info<< " : ";
        if constexpr (UseTypeName)
        {
            printTypeName<base, true>();
        }
        else
        {
            Info<< error::demangle(typeid(base))
                << " (" << sizeof(base) << " bytes)";
        }

        Info<< " cmpt-type=";
        printDataTypeId(UPstream_dataType<Type>::datatype_id);
        Info<< " count=" << count << nl;
    }
}


template<class BinaryOp>
void printOpCodeTraits(BinaryOp bop, std::string_view name)
{
    Info<< "op: " << name << ' ';

    printOpCodeId(UPstream_opType<BinaryOp>::opcode_id);
    Info<< nl;
}


template<class DataType, class BinaryOp>
void printOpCodeTraits(BinaryOp bop, std::string_view name)
{
    Info<< "op: " << name << ' ';

    printOpCodeId(UPstream_opType<BinaryOp>::opcode_id);

    if constexpr (!std::is_void_v<DataType>)
    {
        if constexpr (UPstream_basic_dataType<DataType>::value)
        {
            Info<< " [supported type]";
        }
        else
        {
            Info<< " [disabled]";
        }
    }
    Info<< nl;
}


template<class DataType, class BinaryOp>
void print_data_opType(BinaryOp bop, std::string_view name)
{
    Info<< "op: " << name << ' ';

    printOpCodeId(UPstream_data_opType<BinaryOp, DataType>::opcode_id);

    const bool ok = UPstream_data_opType<BinaryOp, DataType>::value;

    Info<< " okay=" << ok << nl;
}



template<class Type>
int check_simple(std::string_view name = std::string_view())
{
    // The sizing factor is constexpr
    constexpr std::streamsize count = UPstream_dataType<Type>::size(1);

    static_assert
    (
        (count == 1),
        "Code does not (yet) work with aggregate types"
    );

    Info<< "check_simple: " << name << ": " << count << nl;
    return count;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    // contiguous_element must resolve arithmetic types directly
    // (self pass-through), without requiring any pTraits<T>
    // specialization. 'long double' and 'unsigned char' have none.
    static_assert
    (
        std::is_same_v<long double, contiguous_element_t<long double>>
    );
    static_assert
    (
        std::is_same_v<unsigned char, contiguous_element_t<unsigned char>>
    );
    static_assert
    (
        std::is_same_v<bool, contiguous_element_t<bool>>
    );

    printPstreamTraits<bool>();
    printPstreamTraits<label>();

    printPstreamTraits<char, false>("<char>");
    printPstreamTraits<signed char, false>("<signed char>");
    printPstreamTraits<unsigned char, false>("<unsigned char>");

    printPstreamTraits<int8_t, false>("<int8_t>");
    printPstreamTraits<uint8_t, false>("<uint8_t>");

    printPstreamTraits<int16_t, false>("<int16_t>");
    printPstreamTraits<uint16_t, false>("<uint16_t>");

    printPstreamTraits<int>("<int>");
    printPstreamTraits<long>("<long>");
    printPstreamTraits<unsigned>("<unsigned>");
    printPstreamTraits<unsigned int>("<unsigned int>");
    printPstreamTraits<unsigned long>("<long long>");

    printPstreamTraits<const float>();
    printPstreamTraits<floatVector>();

    check_simple<floatVector>("vector<float>");

    printPstreamTraits<Vector2D<float>>();
    printPstreamTraits<Vector2D<double>>();
    printPstreamTraits<boundBox>();
    printPstreamTraits<treeBoundBox>();

    printPstreamTraits<scalar>();
    printPstreamTraits<double>();
    printPstreamTraits<doubleVector>();

    // Avoid typeName for barycentric. It is declared, but not defined
    printPstreamTraits<barycentric, false>("barycentric");

    printPstreamTraits<complex>();     // Uses specialized pTraits_...

    printPstreamTraits<boolVector>();  // Uses specialized pTraits_...
    printPstreamTraits<floatVector>();
    printPstreamTraits<doubleVector>();
    printPstreamTraits<tensor>();
    printPstreamTraits<word>();
    printPstreamTraits<std::string>();

    // This will not identify properly at the moment...
    printPstreamTraits<FixedList<doubleVector, 4>>();
    printPstreamTraits<FixedList<word, 4>>();
    printPstreamTraits<labelPair>();
    printPstreamTraits<Pair<int32_t>>();
    printPstreamTraits<Pair<int64_t>>();
    printPstreamTraits<Pair<float_t>>();
    printPstreamTraits<Pair<double>>();
    printPstreamTraits<Pair<word>>();
    printPstreamTraits<Pair<floatVector>>();
    printPstreamTraits<Pair<doubleVector>>();

    printPstreamTraits<Tuple2<scalar,scalar>>();
    printPstreamTraits<Tuple2<scalar,label>>();
    printPstreamTraits<Tuple2<scalar,vector>>();
    printPstreamTraits<Tuple2<float,floatVector>>();
    printPstreamTraits<Tuple2<floatVector,floatVector>>();
    printPstreamTraits<Tuple2<float,Pair<float>>>();
    printPstreamTraits<Tuple2<float,FixedList<float,3>>>();
    printPstreamTraits<scalarMinMax>();


    Info<< nl
        << "========" << nl
        << "Mapping of binary reduction ops" << nl;

    printOpCodeTraits(minOp<scalar>{}, "min");
    printOpCodeTraits(maxOp<vector>{}, "max");
    printOpCodeTraits(sumOp<vector>{}, "sum");
    printOpCodeTraits(plusOp<vector>{}, "plus");
    printOpCodeTraits(multiplyOp<scalar>{}, "multiply");
    printOpCodeTraits(divideOp<vector>{}, "divide");
    printOpCodeTraits(minMagSqrOp<vector>{}, "minMagSqr");

    printOpCodeTraits(bitAndOp<vector>{}, "bitAnd<vector>");
    printOpCodeTraits(bitOrOp<vector>{}, "bitOr<vector>");

    printOpCodeTraits(bitOrOp<float>{}, "bitOr<float>");
    printOpCodeTraits(bitAndOp<unsigned>{}, "bitAnd<unsigned>");
    printOpCodeTraits(bitOrOp<unsigned>{}, "bitOr<unsigned>");

    printOpCodeTraits<vector>(sumOp<vector>{}, "sum");
    printOpCodeTraits(sumOp<scalarMinMax>{}, "sum");

    printOpCodeTraits(std::plus<>{}, "sum");
    printOpCodeTraits<bool>(std::plus<>{}, "sum");
    printOpCodeTraits<vector>(std::plus<>{}, "sum");

    // Expect success
    Info<< nl << "expect success" << nl;
    print_data_opType<vector>(maxOp<scalar>(), "maxOp(scalar)");
    print_data_opType<unsigned>(bitOrOp<unsigned>(), "bitOrOp(unsigned)");
    print_data_opType<uint8_t>(bitOrOp<uint8_t>(), "bitOrOp(uint8_t)");
    print_data_opType<uint16_t>(bitOrOp<uint16_t>(), "bitOrOp(uint16_t)");

    // Even allow signed integrals
    print_data_opType<int>(bitOrOp<int>(), "bitOrOp(int)");
    print_data_opType<int8_t>(bitOrOp<int8_t>(), "bitOrOp(int8_t)");

    // Failure - supported op, unsupported data type.
    Info<< nl << "expect failure" << nl;
    print_data_opType<bool>(maxOp<scalar>(), "maxOp(scalar, bool)");
    print_data_opType<bool>(bitOrOp<unsigned>(), "bitOrOp(unsigned, bool)");

    // False positives. Failure - supported op, unsupported data type.
    Info<< nl << "false positives" << nl;
    print_data_opType<void>(maxOp<bool>(), "maxOp(bool, void)");
    print_data_opType<float>(bitOrOp<unsigned>(), "bitOrOp(unsigned, float)");

    Info<< nl << "expect failure" << nl;
    print_data_opType<boundBox>(minOp<boundBox>(), "minOp(boundBox)");
    print_data_opType<treeBoundBox>(maxOp<treeBoundBox>(), "maxOp(tbb)");

    print_data_opType<boundBox>(sumOp<boundBox>(), "sumOp(boundBox)");
    print_data_opType<treeBoundBox>(sumOp<treeBoundBox>(), "sumOp(tbb)");


    Info<< nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
