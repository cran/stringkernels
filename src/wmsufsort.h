/* ***** BEGIN LICENSE BLOCK *****
 * Version: MPL 1.1
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * The Original Code is the Suffix Array based String Kernel.
 *
 * The Initial Developer of the Original Code is
 * Statistical Machine Learning Program (SML), National ICT Australia (NICTA).
 * Portions created by the Initial Developer are Copyright (C) 2006
 * the Initial Developer. All Rights Reserved.
 *
 * Contributor(s):
 *
 *   Choon Hui Teo <ChoonHui.Teo@rsise.anu.edu.au>
 *   S V N Vishwanathan <SVN.Vishwanathan@nicta.com.au>
 *
 * ***** END LICENSE BLOCK ***** */


// File    : sask/Code/W_msufsort.h
//
// Authors : Choon Hui Teo      (ChoonHui.Teo@rsise.anu.edu.au)
//           S V N Vishwanathan (SVN.Vishwanathan@nicta.com.au)
//
// Created : 09 Feb 2006
//
// Updated : 24 Apr 2006
//           13 Jul 2007  : use MSufSort v3.1 instead of v2.2

// Wrapper for Michael Maniscalco's MSufSort version 3.1 algorithm
#ifndef W_MSUFSORT_H
#define W_MSUFSORT_H

#include "datatype.h"
#include "isafactory.h"
#include "msufsort.h"


class W_msufsort : public I_SAFactory
{

 public:

	///Variables

	//'Declaration of object POINTERS, no initialization needed.
	//'If Declaration of objects, initialize them in member initialization list.
	MSufSort *msuffixsorter;
	
	///Constructor
	W_msufsort();

	///Destructor
	virtual ~W_msufsort();

	///Methods
	ErrorCode ConstructSA(SYMBOL *text, const UInt32 &len, UInt32 *&array);

};
#endif
