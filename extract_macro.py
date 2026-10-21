#!/usr/bin/env python
"""
Extract a Python macro from a SerialEM settings file into a .py file for testing.

SerialEM stores scripts in the settings file as blocks:

    Macro   <index>
    #!Python
    ...
    EndMacro

<index> is 0-based; the SerialEM GUI shows it as "Script <index+1>".

Usage:
    python extract_macro.py --list
    python extract_macro.py 21                  # by 0-based Macro index
    python extract_macro.py 21 -o test_ronchi.py
    python extract_macro.py -n quick_ronchi     # by ScriptName/MacroName substring
"""

import argparse
import re
import sys

DEFAULT_SETTINGS = 'SerialEMsettings-scripts_AC.txt'
NAME_RE = re.compile(r'^\s*#?\s*(?:ScriptName|MacroName)\s*:?\s*(.*)$')


def parse_macros(path):
    """Return dict of {index: list of lines} for all Macro blocks."""
    macros = {}
    current = None
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if current is None:
                m = re.match(r'^Macro\s+(\d+)\s*$', line)
                if m:
                    current = int(m.group(1))
                    macros[current] = []
            elif line == 'EndMacro':
                current = None
            else:
                macros[current].append(line)
    return macros


def is_python(lines):
    return bool(lines) and lines[0].strip().lower().startswith('#!python')


def macro_name(lines):
    for line in lines[:10]:
        m = NAME_RE.match(line)
        if m and m.group(1).strip():
            return m.group(1).strip()
    return ''


def safe_filename(name, index):
    stem = re.sub(r'[^A-Za-z0-9_]+', '_', name).strip('_') or 'macro'
    return 'macro%d_%s.py' % (index, stem)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('index', nargs='?', type=int,
                        help='0-based Macro index in the settings file (GUI Script number - 1)')
    parser.add_argument('-n', '--name', help='select macro by ScriptName/MacroName substring')
    parser.add_argument('-s', '--settings', default=DEFAULT_SETTINGS,
                        help='SerialEM settings file (default: %(default)s)')
    parser.add_argument('-o', '--output', help='output .py file (default: derived from name)')
    parser.add_argument('-l', '--list', action='store_true', help='list macros and exit')
    parser.add_argument('-f', '--force', action='store_true',
                        help='extract even if the macro is not a Python macro')
    args = parser.parse_args()

    macros = parse_macros(args.settings)

    if args.list or (args.index is None and args.name is None):
        for i in sorted(macros):
            lines = macros[i]
            kind = 'py ' if is_python(lines) else 'sem'
            print('%3d  (Script %2d)  %s  %s' % (i, i + 1, kind, macro_name(lines)))
        return

    if args.index is not None:
        if args.index not in macros:
            sys.exit('Macro %d not found in %s' % (args.index, args.settings))
        index = args.index
    else:
        matches = [i for i in sorted(macros)
                   if args.name.lower() in macro_name(macros[i]).lower()]
        if not matches:
            sys.exit('No macro name matches "%s"' % args.name)
        if len(matches) > 1:
            print('Multiple matches for "%s":' % args.name)
            for i in matches:
                print('%3d  %s' % (i, macro_name(macros[i])))
            sys.exit('Please select by index.')
        index = matches[0]

    lines = macros[index]
    if not is_python(lines) and not args.force:
        sys.exit('Macro %d is not a Python macro (no #!Python line). Use -f to extract anyway.'
                 % index)

    if is_python(lines):
        lines = lines[1:]  # drop the #!Python line
    output = args.output or safe_filename(macro_name(lines) or 'macro', index)
    with open(output, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    print('Wrote macro %d (Script %d) "%s" to %s (%d lines)'
          % (index, index + 1, macro_name(lines), output, len(lines)))


if __name__ == '__main__':
    main()
