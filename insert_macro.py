#!/usr/bin/env python
"""
Paste a Python file back into a macro of a SerialEM settings file.
Counterpart of extract_macro.py.

The body of "Macro <index>" ... "EndMacro" is replaced by "#!Python" followed
by the file contents. If the macro does not exist yet, a new block is added
in index order. A backup of the settings file is written to <settings>.bak.

<index> is 0-based; the SerialEM GUI shows it as "Script <index+1>".
It must be below MaxMacros in the settings file (60 by default); blocks at
MaxMacros and above hold one-line scripts.

SerialEM writes its scripts back to the settings file when it saves settings
or exits, so close SerialEM (or do not save settings) before running this,
then reopen the settings.

Usage:
    python insert_macro.py macro6_quick_ronchi_magic_phase_measurment.py 6
    python insert_macro.py my_script.py 6 -s other_settings.txt
"""

import argparse
import re
import shutil
import sys

DEFAULT_SETTINGS = 'SerialEMsettings-scripts_AC.txt'
DEFAULT_MAX_MACROS = 60
MACRO_RE = re.compile(r'^Macro\s+(\d+)\s*$')


def read_lines(path):
    with open(path, 'rb') as f:
        raw = f.read()
    newline = '\r\n' if b'\r\n' in raw else '\n'
    return raw.decode('utf-8', errors='replace').splitlines(), newline


def max_macros_in(lines):
    """MaxMacros from the settings file, SerialEM default 60 if not found."""
    for line in lines:
        m = re.match(r'^MaxMacros\s+(\d+)\s*$', line)
        if m:
            return int(m.group(1))
    return DEFAULT_MAX_MACROS


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pyfile', help='python file to paste in')
    parser.add_argument('index', type=int,
                        help='0-based Macro index in the settings file (GUI Script number - 1)')
    parser.add_argument('-s', '--settings', default=DEFAULT_SETTINGS,
                        help='SerialEM settings file (default: %(default)s)')
    parser.add_argument('--no-backup', action='store_true', help='do not write <settings>.bak')
    args = parser.parse_args()

    lines, newline = read_lines(args.settings)
    max_macros = max_macros_in(lines)
    if not 0 <= args.index < max_macros:
        sys.exit('Macro index must be 0 to %d (MaxMacros %d in %s). Higher blocks are one-line scripts.'
                 % (max_macros - 1, max_macros, args.settings))
    body, _ = read_lines(args.pyfile)
    if any(line == 'EndMacro' for line in body):
        sys.exit('%s contains a line "EndMacro", which would end the macro early' % args.pyfile)
    if not (body and body[0].strip().lower().startswith('#!python')):
        body = ['#!Python'] + body

    # locate existing block and the insertion point for a new one
    start = end = None
    insert_at = None
    last_end = None
    current = None
    for i, line in enumerate(lines):
        m = MACRO_RE.match(line)
        if current is None and m:
            current = int(m.group(1))
            if current == args.index:
                start = i
            elif current > args.index and insert_at is None:
                insert_at = i
        elif current is not None and line == 'EndMacro':
            if current == args.index:
                end = i
            last_end = i
            current = None

    if start is not None:
        if end is None:
            sys.exit('Macro %d has no EndMacro in %s' % (args.index, args.settings))
        old_len = end - start - 1
        new_lines = lines[:start + 1] + body + lines[end:]
        action = 'Replaced macro %d (%d lines -> %d lines)' % (args.index, old_len, len(body))
    else:
        if insert_at is None:
            if last_end is None:
                sys.exit('No Macro blocks found in %s' % args.settings)
            insert_at = last_end + 1
        block = ['Macro\t%d' % args.index] + body + ['EndMacro']
        new_lines = lines[:insert_at] + block + lines[insert_at:]
        action = 'Added new macro %d (%d lines)' % (args.index, len(body))

    if not args.no_backup:
        shutil.copy2(args.settings, args.settings + '.bak')
    with open(args.settings, 'w', encoding='utf-8', newline='') as f:
        f.write(newline.join(new_lines) + newline)
    print('%s (Script %d) in %s from %s' % (action, args.index + 1, args.settings, args.pyfile))


if __name__ == '__main__':
    main()
