import click
import re
import os
from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
from . import __version__
from .core import Analyze


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.group(context_settings=CONTEXT_SETTINGS, chain=True)
@click.pass_context
@optgroup.group('general options')
@optgroup.option('-g', '--genecode', 'genecode', metavar='<int>', default=1, type=click.INT,
                 show_default=True, help='This option is use for change codon table.')
def cli(ctx, genecode):
    ctx.ensure_object(dict)
    ctx.obj['code'] = genecode
    pass


@cli.command('analyze',help=r'Calculate indexes.Try cubat analyze -h for more...')
@click.pass_context
@click.argument('input', type=click.Path(exists=True), metavar='INPUT')
@optgroup.group('CUB analyze', help='computing indexes for gene.')
@optgroup.option('-e', '--enc', 'enc', is_flag=True, help=r'a flag option to compute ite.')
@optgroup.option('-x', '--chi^2', 'chi2', is_flag=True, help=r'a flag option to compute chi-squared.')
@optgroup.option('-c', '--cai', 'cai', is_flag=True, help=r'a flag option to compute cai.')
@optgroup.option('--cr', '--cairef', 'cairef', type=click.Path(),
                 metavar='<reference>', help=r'reference for cai computing, refering to cubat\example\cai_ref.csv.')
@optgroup.option('-f', '--fop', 'fop', is_flag=True, help='a flag option to compute fop.')
@optgroup.option('--fo', '--fopopt', 'fopopt', type=click.Path(),
                 metavar='<optimal>',
                 help=r'optimal codons for fop computing, refering to cubat\example\fop_opt.csv.')
@optgroup.option('-b', '--cbi', 'cbi', is_flag=True, help='a flag option to compute cbi.')
@optgroup.option('--bo', '--cbiopt', 'cbiopt', type=click.Path(),
                 metavar='<optimal>',
                 help=r'optimal codons for cbi computing, refering to cubat\example\cbi_opt.csv.')
@optgroup.option('-i', '--ite', 'ite', is_flag=True, help='a flag option to compute ite.')
@optgroup.option('--ir', '--iteref', 'iteref', type=click.Path(),
                 metavar='<reference>', help=r'reference for ite computing, refering to cubat\example\ite_ref.csv.')
@optgroup.option('-t', '--tai', 'tai', is_flag=True, help='a flag option to compute tai.')
@optgroup.option('--tg', '--trnagcn', 'taigcn', type=click.Path(),
                 metavar='<trna copy>',
                 help=r'trna gene copy number for tai computing, refering to cubat\example\tai_gcn.csv.')
@optgroup.option('--ts', '--tais', 'tais', type=click.Path(),
                 metavar='<constraint>',
                 help=r'selective constraint for tai computing, refering to cubat\example\tai_s.csv.'
                      r'if none, the default value for eucaryon will be used')


# @optgroup.group('codon analyze', help='computing indexes for codon.')
@optgroup.option('-r', '--rscu', 'rscu', is_flag=True, help=r'a flag option to compute rscu.')
@optgroup.option('-s', '--csc', 'csc', is_flag=True, help=r'a flag option to compute csc.')
@optgroup.option('--ml', '--mrnahl', 'mrnahl', type=click.Path(),
                 metavar='<mrna half-life>',
                 help=r'half-life of mrna for csc computing, refering to cubat\example\csc_mranhl.csv.')

@optgroup.group('output and saving')
@optgroup.option('-o', '--output', 'output', help=r'output path.')
@optgroup.option('-p', '--prefix', 'prefix', default='', help=r'prefix for output files.')
@optgroup.option('-a', '--save', is_flag=True, help=r'a flag option of saving.')
def cubat_analyze(ctx, enc, chi2, cai, cairef, fop, fopopt, cbi, cbiopt, ite, iteref, tai, taigcn, tais, rscu, csc,
                  mrnahl, output, save, prefix, input):
    # process input file. get format.
    ctx.obj['input'] = input
    dire = False
    if os.path.isdir(input):
        dire = True
        format_list = []
        file_list = []
        for filepath, dirnames, filenames in os.walk(input):
            for filename in filenames:
                if re.findall('\w{1,}$', filename)[0] in ['fasta']:
                    # todo 注意其它文件格式
                    format_list.append(re.findall('\w{1,}$', filename)[0])
                    file_list.append(os.path.join(filepath, filename))

    else:
        # todo 采用字典解决
        if re.findall('\w{1,}$', input)[0] == 'FAS':
            format_list = ['fasta']
        else:
            format_list = re.findall('\w{1,}$', input)
        file_list = [input]

    ctx.obj['file_list']=file_list

    # process other parameters
    try:
        genecode = int(ctx.obj['code'])
    except:
        raise AttributeError('please give an integer for genecode')

    # do the computing
    for amount in range(len(file_list)):
        result = Analyze(file_list[amount], file_format=format_list[amount], genecode=genecode, enc=enc, cai=cai,
                         cai_ref=cairef, X2=chi2, fop=fop, fop_opt=fopopt,
                         cbi=cbi, cbi_opt=cbiopt, tai_s=tais,
                         tai=tai, tai_gcn=taigcn, rscu=rscu, csc=csc, csc_ref=mrnahl, output=output, save=save, prefix=prefix)
        ctx.obj[file_list[amount]] = result
    #todo 生成器


    pass


@cli.command('plot', help=r'Try cubat plot -h for more...')
@click.pass_context
@click.argument('input', type=click.Path(), metavar='INPUT', default='')
@optgroup.group('plots')
@optgroup.option('--rb', '--rscubar', 'rscubar', is_flag=True, help=r'a flag option to output a barplot for rscu.')
@optgroup.group('output')
@optgroup.option('-o', '--output', 'output', help=r'output path.')
@optgroup.option('-p', '--prefix', 'prefix', default='', help=r'prefix for output files.')
def cubat_plot(ctx, rscubar, output, prefix, input):
    genecode = ctx.obj['code']

    rscu = False
    if rscubar:
        rscu = True

    if ctx.obj['input'] == input or input == '':
        for files in ctx.obj['file_list']:
            ctx.obj[files].plot(rscu_barplot=rscubar, output=output, prefix=prefix)
    else:
        cubat_analyze(input=input, genecode=genecode, output=output, prefix=prefix, rscu=rscu)
        for files in ctx.obj['file_list']:
            ctx.obj[files].plot(rscu_barplot=rscubar, output=output, prefix=prefix)
    return

if __name__ == '__main__':
    cli()
