//
// Created by krab1k on 2019-11-14.
//

#pragma once

#include <string>

std::string to_lowercase(const std::string &from);

std::string to_uppercase(const std::string &from);

bool starts_with(const std::string &text, const std::string &prefix);

bool ends_with(std::string const &fullString, std::string const &suffix);

std::string trim(const std::string &text);
