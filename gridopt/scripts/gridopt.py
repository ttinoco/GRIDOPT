#*****************************************************#
# This file is part of GRIDOPT.                       #
#                                                     #
# Copyright (c) 2015, Tomas Tinoco De Rubira.         #
#                                                     #
# GRIDOPT is released under the BSD 2-clause license. #
#*****************************************************#

from __future__ import print_function
import sys
import pstats
import argparse
import cProfile
from ..power_flow import new_method, PFmethodError

methods = ['ACOPF', 'ACPF', 'DCOPF', 'DCPF']

def create_parser():
    
    # Construct parser
    parser = argparse.ArgumentParser(description='Power grid optimization package.')
    
    # Case
    parser.add_argument('case', type=str,
                        help='filename of power flow case to solve')

    # Method
    parser.add_argument('method', choices=methods,metavar='method',
                        help='PF or OPF method')
    
    # Params
    class ParamAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            for value in values:
                if '=' not in value:
                    raise argparse.ArgumentTypeError("invalid parameter name-value pair")
                n,v = value.split('=')
                params = getattr(namespace, self.dest)
                setattr(namespace, self.dest,params+[(n, v)])
    parser.add_argument('--parameters', action=ParamAction, dest='params', nargs='*', default=[],
                        help='parameter name-value pairs')

    # Profile
    parser.add_argument('--profile', action='store_true', default=False,
                        help='flag for profiling execution')

    # Flat start
    parser.add_argument('--flatstart', action='store_true', default=False,
                        help='flag for flat starting point')

    # Quiet
    parser.add_argument('--quiet', action='store_true', default=False,
                        help='flag for supressing output')
    
    # Write
    parser.add_argument('--write',
                        dest='outfile',
                        type=str,
                        default='',
                        help='name of output file')

    return parser

def main(args=None):

    import pfnet

    parser = create_parser()
    args = parser.parse_args(args=args if args is not None else sys.argv[1:])
    
    try:
        
        # Network
        net = pfnet.Parser(args.case).parse(args.case)
        if not args.quiet:
            net.show_components()
        
        if args.flatstart:

            # Flat start
            for bus in net.buses:
                bus.v_mag = 1
                bus.v_ang = 0
        
        # Method
        method = new_method(args.method)
            
        # Parameters
        method.set_parameters(strparams=dict(args.params))
        if args.quiet:
            method.set_parameters({'quiet': True})
            
        if args.profile:

            # Profile
            cProfile.runctx('method.solve(net)', globals(), locals(), '.prof')
            pstats.Stats('.prof').strip_dirs().sort_stats('cumulative').print_stats(20)
        else:

            # Solve
            try:
                method.solve(net)
                if not args.quiet:
                    print(method.results['solver status'])
            except PFmethodError as e:
                if not args.quiet:
                    print(e)
            
            # Update
            method.update_network(net)
                
            # Show network props
            if not args.quiet:
                net.show_properties()

            # Write
            if args.outfile:
                p = pfnet.Parser(args.outfile)
                p.write(net, args.outfile)

    finally:
        pass
    
# Main function            
if __name__ == "__main__":
    main()
