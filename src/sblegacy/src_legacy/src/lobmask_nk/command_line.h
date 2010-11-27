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

/* ///////////////////////////////////////////////////////////////////
//
//  command_line.h
//    Functions and classes for handline Command Line input
//
//////////////////////////////////////////////////////////////////////
//
//  Gal Sela
//  --------------------------------------------------
//  $Id: command_line.h,v 1.6 1998/11/19 18:26:14 sela Exp $
//
/////////////////////////////////////////////////////////////////// */

#ifndef GS_COMMANDLINE_H
#define GS_COMMANDLINE_H

#include <string.h>

#ifdef __cplusplus
namespace gs {
#endif

/* -----------------------------------------------------------
//  commandLine - extracts options from command line
//
//    o Looks for tag in argv - returns argv[tag_pos+offset]
//    o tag can be a partial match to argv[] string
//    o returns NULL if no match found
*/
#ifdef __cplusplus
char *commandLine(int argc, char **argv, char *tag, int offset = 0)
#else
char *commandLine(int argc, char **argv, char *tag, int offset)
#endif
{
    int j;

    for(j=0; j<argc; j++){
	if(!strncmp(argv[j], tag, strlen(tag))){
	    if(j+offset < argc){
		return(argv[j+offset]);
	    }
	    else{
		return(NULL);
	    }
	}
    }
    
    return(NULL);
}

#ifdef __cplusplus
/////////////////////////////////////////////////////////////////
//
//  CommandLine class
//    - encapsulates CommandLine look-up
//
class CommandLine
{
  public:
    CommandLine(int argc, char **argv)
      : _argc(argc), _argv(argv), _current(0) { }

    CommandLine &find(const char *tag, int offset = 0) {
      _current = _argc;
      for(int j=0; j<_argc; j++){
	if(!strncmp(_argv[j], tag, strlen(tag))){
	  if(j+offset < _argc) _current = j+offset;
	  break;
	}
      }
      return(*this);
    }
    
    operator int() { return(_current < _argc); }
    operator const char *() { return(_argv[_current]); }

    CommandLine &operator++() {
      ++_current;
      if(_current > _argc) _current = _argc;
      return(*this);
    }
  
    CommandLine operator++(int) {
      CommandLine tmp = *this;
      ++*this;
      return(tmp);
    }
 
    char * operator*() {
      return(_argv[_current]);
    }

    char * operator[](int offset) {
      return(_argv[_current+offset]);
    }    

    int current() const { return(_current); }

  private:
    int _argc;
    char **_argv;
    int _current;
};

#endif // __cplusplus

#ifdef __cplusplus
} // namespace gs
#endif


#endif // GS_COMMANDLINE_H


/* ///////////////////////////////////////////////////////////////////
//
//  $Log: command_line.h,v $
//  Revision 1.6  1998/11/19 18:26:14  sela
//  o Add current() to class CommandLine
//
//  Revision 1.5  1998-10-28 16:19:23-05  sela
//  o Make default offset = 0 in find()
//  o Add automatic conversion to const char *
//
//  Revision 1.4  1998-10-28 16:05:33-05  sela
//  o Move string.h before namespace gs
//
//  Revision 1.3  1998-10-08 16:31:33-04  sela
//  o convert to stdc++
//
//  Revision 1.2  1998-08-19 14:47:46-04  sela
//  o Fix bug in setting offset of find() command
//  o Set default value for offset in find() to 1
//
//  Revision 1.1  1998-08-19 11:25:05-04  sela
//  Initial revision
//
//
/////////////////////////////////////////////////////////////////// */




