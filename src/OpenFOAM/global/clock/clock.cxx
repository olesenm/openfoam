/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "clock.H"

#include <cstdio>
#include <ctime>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

// static inline const std::tm* get_gmtime()
// {
//     std::time_t t = std::time(nullptr);  // ie, now()
//     return std::gmtime(&t);
// }


static inline const std::tm* get_localtime()
{
    std::time_t t = std::time(nullptr);  // ie, now()
    return std::localtime(&t);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

std::string Foam::clock::dateTime()
{
    // ISO-8601 format (yyyy-mm-ddThh:mm:ss)
    constexpr auto fmt_string = "%04d-%02d-%02dT%02d:%02d:%02d";
    constexpr int  fmt_width  = 19;
    std::string str(fmt_width, '\0');

    const auto* curr = get_localtime();

    // if (!curr)
    // {
    //     return "1970-01-01T00:00:00";
    // }

    std::snprintf
    (
        str.data(),
        str.size()+1,
        fmt_string,
        (curr->tm_year + 1900), (curr->tm_mon + 1), (curr->tm_mday),
        curr->tm_hour, curr->tm_min, curr->tm_sec
    );

    return str;
}


std::string Foam::clock::iso_date()
{
    // ISO-8601 format (yyyy-mm-dd)
    constexpr auto fmt_string = "%04d-%02d-%02d";
    constexpr int  fmt_width  = 10;
    std::string str(fmt_width, '\0');

    const auto* curr = get_localtime();

    // if (!curr)
    // {
    //     return "1970-01-01";
    // }

    std::snprintf
    (
        str.data(),
        str.size()+1,
        fmt_string,
        (curr->tm_year + 1900), (curr->tm_mon + 1), (curr->tm_mday)
    );

    return str;
}


std::string Foam::clock::date()
{
    // Date format (MON dd yyyy)
    constexpr auto fmt_string = "%.3s %02d %04d";
    constexpr int  fmt_width  = 11;
    std::string str(fmt_width, '\0');

    const auto* curr = get_localtime();

    // if (!curr)
    // {
    //     return "Jan 01 1970";
    // }

    constexpr char months[] =
    (
        "Jan" "Feb" "Mar" "Apr" "May" "Jun"
        "Jul" "Aug" "Sep" "Oct" "Nov" "Dec"
    );

    std::snprintf
    (
        str.data(),
        str.size()+1,
        fmt_string,
        (months + 3*(curr->tm_mon)), (curr->tm_mday), (curr->tm_year + 1900)
    );

    return str;
}


std::string Foam::clock::clockTime()
{
    // Time format (hh:mm:ss)
    constexpr auto fmt_string = "%02d:%02d:%02d";
    constexpr int  fmt_width  = 8;
    std::string str(fmt_width, '\0');

    const auto* curr = get_localtime();

    // if (!curr)
    // {
    //     return "00:00:00";
    // }

    std::snprintf
    (
        str.data(),
        str.size()+1,
        fmt_string,
        curr->tm_hour, curr->tm_min, curr->tm_sec
    );

    return str;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clock::clock()
:
    start_(now()),
    last_(start_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

double Foam::clock::elapsedClockTime() const
{
    last_ = now();
    return std::difftime(last_, start_);
}


double Foam::clock::clockTimeIncrement() const
{
    const auto prev(last_);
    last_ = now();
    return std::difftime(last_, prev);
}


// ************************************************************************* //
