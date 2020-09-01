/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 - 2020 Iowa State University
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
//////////////////////////////////////////////////////////////////////////////////

#ifndef PRS_CENCODE_H
#define PRS_CENCODE_H

/***
cencode.h - c source to a base64 encoding algorithm implementation

This is part of the libb64 project, and has been placed in the public domain.
For details, see <a href = "http://sourceforge.net/projects/libb64"> http://sourceforge.net/projects/libb64</a>
*/


typedef enum
{
  step_A, step_B, step_C
} 	base64_encodestep;

typedef struct
{
  base64_encodestep step;
  char result;
  int stepcount;
} base64_encodestate;

inline void base64_init_encodestate(base64_encodestate* state_in);

inline char base64_encode_value(char value_in);

inline int base64_encode_block(const char* plaintext_in, int length_in, char* code_out, base64_encodestate* state_in);

inline int base64_encode_blockend(char* code_out, base64_encodestate* state_in);



const int CHARS_PER_LINE = 72;

inline void base64_init_encodestate(base64_encodestate* state_in)
{
  state_in->step = step_A;
  state_in->result = 0;
  state_in->stepcount = 0;
}

inline char base64_encode_value(char value_in)
{
  static const char* encoding = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  if (value_in > 63) return '=';
  return encoding[(int)value_in];
}

inline int base64_encode_block(const char* plaintext_in, int length_in, char* code_out, base64_encodestate* state_in)
{
  const char* plainchar = plaintext_in;
  const char* const plaintextend = plaintext_in + length_in;
  char* codechar = code_out;
  char result;
  char fragment;

  result = state_in->result;

  switch (state_in->step)
  {
    while (1)
    {
      case step_A:
        if (plainchar == plaintextend)
        {
          state_in->result = result;
          state_in->step = step_A;
          return codechar - code_out;
        }
      fragment = *plainchar++;
      result = (fragment & 0x0fc) >> 2;
      *codechar++ = base64_encode_value(result);
      result = (fragment & 0x003) << 4;
      case step_B:
        if (plainchar == plaintextend)
        {
          state_in->result = result;
          state_in->step = step_B;
          return codechar - code_out;
        }
      fragment = *plainchar++;
      result |= (fragment & 0x0f0) >> 4;
      *codechar++ = base64_encode_value(result);
      result = (fragment & 0x00f) << 2;
      case step_C:
        if (plainchar == plaintextend)
        {
          state_in->result = result;
          state_in->step = step_C;
          return codechar - code_out;
        }
      fragment = *plainchar++;
      result |= (fragment & 0x0c0) >> 6;
      *codechar++ = base64_encode_value(result);
      result  = (fragment & 0x03f) >> 0;
      *codechar++ = base64_encode_value(result);

      ++(state_in->stepcount);
    }
  }
  /* control should not reach here */
  return codechar - code_out;
}

inline int base64_encode_blockend(char* code_out, base64_encodestate* state_in)
{
  char* codechar = code_out;

  switch (state_in->step)
  {
    case step_B:
      *codechar++ = base64_encode_value(state_in->result);
      *codechar++ = '=';
      *codechar++ = '=';
      break;
    case step_C:
      *codechar++ = base64_encode_value(state_in->result);
      *codechar++ = '=';
      break;
    case step_A:
      break;
  }

  return codechar - code_out;
}

#endif //PRS_CENCODE_H
