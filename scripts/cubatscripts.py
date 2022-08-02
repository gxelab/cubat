import warnings
import click
import os
import sys
sys.path.append('C:/Users/YuanYe/PycharmProjects/pythonProject/projects/venv/Lib')
root_path = os.path.abspath(os.path.dirname(__file__))
main_path=os.path.abspath(os.path.join(root_path,os.path.pardir))
sys.path.append(os.path.abspath(main_path))
import main



@click.group(chain=True)
@click.option('-i', '--input', 'input', type=click.Path(),
              help='Files you want to analyze. Welcomes fasta files or a dictionary.'
    ,multiple=True, default=())
@click.option('-p', '--parameter', 'parameter', type=click.Path(), default='')
@click.option('-o', '--output', 'output', type=click.Path(), default='')

@click.pass_context
def cli(ctx,input, parameter, output):
    # todo exists=True
    # judge whether the input file is a dir
    ctx.ensure_object(dict)
    ctx.obj['selected']=False
    ctx.obj['dire'] = False
    infile = input
    if os.path.isdir(infile[0]):
        ctx.obj['dire'] = True
        if len(infile)>1:
            warnings.warn('Only the first directory would be read. You may merge them into one')
    ctx.obj['inp'] = infile

    # judge whether parameter file and putput path are existed.
    parafile = parameter
    if parafile != '':
        if not os.path.exists(parafile):
            raise ValueError('the parameter file did not exist.')
        else:
            ctx.obj['para'] = parafile
    else:
        ctx.obj['para'] = parafile

    outpath = output
    if outpath != '':
        if not os.path.exists(outpath):
            raise ValueError('the output path did not exist.')
        else:
            ctx.obj['out'] = outpath
    else:
        ctx.obj['out'] = outpath

    pass


@cli.command('1', help='select_parameter')
@click.pass_context
def select_parameters(ctx):
    # ctx.ensure_object(dict)
    sel_class = main.select_parameters()
    ctx.obj['par_class'] = sel_class
    ctx.obj['selected'] = True

    return



@cli.command('2', help='codon_bias_analyze')
@click.option('genecode', '--genecode','-c', default=0, type=int)
@click.pass_context
def codon_bias_analyze(ctx,genecode):
    files_list = []
    # todo 注意部分文件相对部分绝对的情况
    # process file paths passed by preceding part.
    # process the input
    if ctx.obj['dire']:
        if os.path.isabs(ctx.obj['inp'][0]):
            inpath = ctx.obj['inp'][0]
        else:
            inpath = os.getcwd() + '/' + ctx.obj['inp'][0]
        for filepath, dirnames, filenames in os.walk(inpath):
            for filename in filenames:
                if 'fasta'or'FAS' in filename:
                    files_list.append(filepath + '/' + filename)
                if 'parameter_file.xlsx' in filename:
                    ctx.obj['para'] = filepath + '/' + 'parameter_file.xlsx'
    else:
        if os.path.isabs(ctx.obj['inp'][0]):
            inpath = ctx.obj['inp']
        else:
            inpath = []
            for path in ctx.obj['inp']:
                inpath.append(os.getcwd() + '/' + path)
        files_list = inpath

    # process genetic code
    if genecode==0:
        while True:
            try:
                genecode=int(click.prompt('please enter the genetic code'))
            except ValueError:
                print('please enter a number')
            if genecode:
                break


    # process the parameter file
    if ctx.obj['selected']:
        if ctx.obj['para']:
            ctx.obj['par_class'].process_parameter_file(ctx.obj['para'])
            parameters = ctx.obj['par_class']
        else:
            parameters = ctx.obj['par_class']

    else:
        if ctx.obj['para']:
            parameters = main.select_parameters(parameter_file=ctx.obj['para'], skip_select=True)

        else:
            if click.confirm("Since you neither selected parameters nor enter the parameter file,"
                             "\ndo you want to select parameters now"):
                parameters = main.select_parameters()
            else:
                parameters = main.select_parameters(skip_select=True)

    #choose indexes you want to analyze.
    try:
        indexes = click.prompt("Press the serial option of indexes you want to compute: "
                                             "\n(1) enc "
                                             "\n(2) ite "
                                             "\n(3) chi_squared "
                                             "\n(4) tai "
                                             "\n(5) fop "
                                             "\n(6) cbi"
                                             "\n(7) cai"
                                             "\n(8) rscu"
                                             "\n(9) csc"
                                             "\n(0) all_indexes(less csc for its unique request for half_life of mRNA)"
                                             "\n(q) quality control"
                                             "\nyour choices").split()
    except ValueError:
        #todo 注意这里的报错。
        print('"Please enter a number"')

    flag_dict = {'1': False, '2': False, '3': False, '4': False, '5': False, '6': False, '7': False, '8': False,
                 '9': False, '0': False, 'q': False}
    for option in indexes:
        flag_dict[option] = True
        if '0' in indexes:
            flag_dict['q'] = True
            for number in range(9):
                flag_dict[str(number)] = True

    #compute indexes for files one by one
    for file in files_list:
        main.codon_bias_analyze(file, parameters=parameters, genecode=genecode, quality_control=flag_dict['q'],
                                enc=flag_dict['1'],
                                ite=flag_dict['2'], X2=flag_dict['3'], tai=flag_dict['4'],
                                fop=flag_dict['5'], cbi=flag_dict['6'], cai=flag_dict['7'], rscu=flag_dict['8'],
                                csc=flag_dict['9'])
    return

#todo 1.特殊参数的修改 2.完善报错与警告 3.输出文件的问题 4.路径 5.参数运算与运算报错 6.注意部分非标准密码子表的问题，参考dambe。 7.可能的话边输入命令行选择边计算？
#todo 8.注意导入文件格式 9.注意命令行与解释器功能不一致的问题


if __name__ == '__main__':
    cli()
