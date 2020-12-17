//
// Created by roryh on 27/10/2020.
//

#include <iostream>
#include <EpiGraph/EpiGraph.hpp>

namespace EpiGraph {
    auto print_banner() -> void {
        std::cout <<
                  " _____       _  ____                 _    \n" <<
                  "| ____|_ __ (_)/ ___|_ __ __ _ _ __ | |__ \n" <<
                  "|  _| | '_ \\| | |  _| '__/ _` | '_ \\| '_ \\\n" <<
                  "| |___| |_) | | |_| | | | (_| | |_) | | | |\n" <<
                  "|_____| .__/|_|\\____|_|  \\__,_| .__/|_| |_|\n" <<
                  "      |_|                     |_|          \n";
    }

    auto print_banner2() -> void {
        std::cout << "                                                                               \n"
                     "                       ,,                                        ,,            \n"
                     "`7MM\"\"\"YMM             db   .g8\"\"\"bgd                          `7MM            \n"
                     "  MM    `7                .dP'     `M                            MM            \n"
                     "  MM   d   `7MMpdMAo.`7MM dM'       ` `7Mb,od8 ,6\"Yb. `7MMpdMAo. MMpMMMb.      \n"
                     "  MMmmMM     MM   `Wb  MM MM            MM' \"'8)   MM   MM   `Wb MM    MM      \n"
                     "  MM   Y  ,  MM    M8  MM MM.    `7MMF' MM     ,pm9MM   MM    M8 MM    MM      \n"
                     "  MM     ,M  MM   ,AP  MM `Mb.     MM   MM    8M   MM   MM   ,AP MM    MM      \n"
                     ".JMMmmmmMMM  MMbmmd' .JMML. `\"bmmmdPY .JMML.  `Moo9^Yo. MMbmmd'.JMML  JMML.    \n"
                     "             MM                                         MM                     \n"
                     "           .JMML.                                     .JMML.                   \n"
                     "";
    }
}