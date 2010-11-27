/*=========================================================================
  
  Copyright 2002, 2006-2010, Dr. Sandra Black
  Linda C. Campbell Cognitive Neurology Unit
  Sunnybrook Health Sciences Center
  
  This file is part of the Sunnybrook Image Software Processing (SIPS) package

  SIPS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

//This Class was written by : Dr. Azhar Quddus
// Nov 20, 2002.
// All rights reserved. Dept. of Cognitive Neurology, Sunnybrook & Womens'


struct list 
{
        unsigned short tail, next, root;
};
  
/*
        Merge two lists with different roots into a single list

        Example:
                Merges two lists with root 1  and root 5.  Lists with root 2
                and 9 remain unchanged.
                                                 before                            after
Address                       Tail    pointer Root                    Tail    pointer Root
        1                       6             4             1          8       4             1
        2                       3             3             2          3       3             2
        3                       3             0             2          3       0             2
        4                       4             6             1          4       6             1
        5                       8             7             5          8       7             1
        6                       6             0             1          6       5             1
        7                       7             8             5          7       8             1
        8                       8             0             5          8       0             1
        9                       9             0             9          9       0             9

Notes:
        The root is the address of the start of the list.

        Tail is pointer to the tail of the list. Only the root node tail
        pointer is valid.

        Also pointer is a pointer to the address of the next element in
        the list. A pointer of 0 means no next element.

        No pointer should       point to an element with a different root
        (this is what separates the multiple lists).

        entry:
             labels points to array of lists which contains
             the lists are in the form;

                 labels[address].tail = pointer to address of last element
                                        in the list. This pointer is only
                                        for the root of the list.
                 labels[address].next = pointer to address of next element
                                        in the list.
                 labels[address].root = root address for list.
                 label_first = root of first list
                 label_last = root of second list

        Lists are merged on exit.
*/
void link(struct list *labels, unsigned short label_first,
        unsigned short label_last)
{
        unsigned short label_current, root;

        if( label_first > label_last )  // Sort the two lists into a higher and
        {                               // lower root
                root = label_last;      
                label_last = label_first;
        }
        else
                root = label_first;

        // walk the higher list to change the root link to new root
        label_current = label_last;

        do
        {
                labels[label_current].root = root;
        } while( 0 != (label_current = labels[label_current].next) );

        //  Connect higher list with lower list
        labels[labels[root].tail].next = label_last;
        labels[root].tail = labels[label_last].tail;
}
