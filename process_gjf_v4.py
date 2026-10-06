#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import sys
from pathlib import Path


DEFAULT_HEADER = """%nprocshared={nproc}
%chk={chk}
%mem={mem}
#p {opt_keyword} {opt_basis_keyword} {functional} {dispersion_keyword}

{title}
"""

DEFAULT_SCF_XQC = """%nprocshared={nproc}
%chk={scf_chk}
%mem={mem}
#p {opt_keyword} {opt_basis_keyword} {functional} {dispersion_keyword} scf=xqc

Initial geometry optimization with SCF=XQC
"""

DEFAULT_OPT_LINK1 = """--Link1--
%OldChk={scf_chk}
%chk={chk}
%nprocshared={nproc}
%mem={mem}
#p {opt_basis_keyword} {functional} {dispersion_keyword} geom=check guess=read {opt_keyword}

Geometry optimization after SCF=XQC

{charge_mult}

{opt_basis_block}
"""

DEFAULT_BASIS = """Fe O 0
{genecp_feo_basis}
****
C H N 0
{genecp_chn_basis}
****

Fe O 0
{genecp_feo_ecp}
"""

DEFAULT_STABLE = """--Link1--
%OldChk={chk}
%chk={stable_chk}
%nprocshared={nproc}
%mem={mem}
#p {stable_basis_keyword} {functional} {dispersion_keyword} geom=check guess=read stable=opt

Wavefunction stability optimization

{charge_mult}

{stable_basis_block}
"""

DEFAULT_FREQ = """--Link1--
%OldChk={post_opt_chk}
%chk={freq_chk}
%nprocshared={nproc}
%mem={mem}
#p freq {freq_basis_keyword} {functional} {dispersion_keyword} geom=check guess=read

Frequency calculation

{charge_mult}

{freq_basis_block}
"""

DEFAULT_SOL = """--Link1--
%OldChk={post_opt_chk}
%chk={sol_chk}
%nprocshared={nproc}
%mem={mem}
#p scrf=(smd,solvent={solvent}) {sol_basis} {functional} {dispersion_keyword} geom=check guess=read {sol_stable_keyword}

SMD single-point calculation

{charge_mult}
"""

# Standalone jobs used by the pure Freq / pure SMD / Freq→SMD modes.
DEFAULT_FREQ_ONLY = """%nprocshared={nproc}
%chk={freq_chk}
%mem={mem}
#p freq {freq_basis_keyword} {functional} {dispersion_keyword}

Frequency calculation
"""

DEFAULT_SOL_ONLY = """%nprocshared={nproc}
%chk={sol_chk}
%mem={mem}
#p scrf=(smd,solvent={solvent}) {sol_basis} {functional} {dispersion_keyword} {sol_stable_keyword}

SMD single-point calculation
"""

DEFAULT_SOL_FROM_FREQ = """--Link1--
%OldChk={freq_chk}
%chk={sol_chk}
%nprocshared={nproc}
%mem={mem}
#p scrf=(smd,solvent={solvent}) {sol_basis} {functional} {dispersion_keyword} geom=check guess=read {sol_stable_keyword}

SMD single-point calculation after Freq

{charge_mult}
"""


def read_template(path_value, default):
    """Only read an external template when the user explicitly provides one."""
    if path_value is None:
        return default

    path = Path(path_value).expanduser()
    if not path.exists():
        raise FileNotFoundError("模板文件不存在：{}".format(path))

    return path.read_text(encoding="utf-8")


def split_sections(lines):
    """
    Split a Gaussian input into:
      route/link0, title, charge+mult+coordinates, tail.

    Blank lines immediately after the charge/multiplicity line are tolerated.
    Coordinates are still read and retained instead of being mistaken for tail.
    """

    n_lines = len(lines)

    def find_blank(start):
        for index in range(start, n_lines):
            if lines[index].strip() == "":
                return index
        return -1

    first_blank = find_blank(0)
    if first_blank == -1:
        raise ValueError("未找到 route section 后的空行。")

    title_start = first_blank + 1
    second_blank = find_blank(title_start)
    if second_blank == -1:
        raise ValueError("未找到标题后的空行。")

    title = lines[title_start:second_blank]

    charge_index = second_blank + 1
    while charge_index < n_lines and lines[charge_index].strip() == "":
        charge_index += 1

    if charge_index >= n_lines:
        raise ValueError("未找到电荷和多重度。")

    charge_mult = lines[charge_index].strip()
    if not re.match(r"^[+-]?\d+\s+[+-]?\d+(?:\s+.*)?$", charge_mult):
        raise ValueError("无法识别电荷/多重度行：{}".format(lines[charge_index]))

    coord_start = charge_index + 1
    while coord_start < n_lines and lines[coord_start].strip() == "":
        coord_start += 1

    # 有些旧输入文件会意外连续写两次相同的电荷/多重度，例如：
    #
    #   1 2
    #   1 2
    #   Fe ...
    #
    # 第一行已经被识别为 charge_mult；这里自动跳过后面紧邻且完全相同
    # 的重复行，避免把第二个 "1 2" 误当成坐标写回新输入。
    while (
        coord_start < n_lines
        and lines[coord_start].strip() == charge_mult
    ):
        coord_start += 1

        # 同时容忍重复行之间夹有空行。
        while coord_start < n_lines and lines[coord_start].strip() == "":
            coord_start += 1

    if coord_start >= n_lines:
        raise ValueError("电荷/多重度后未找到坐标。")

    coord_end = coord_start
    while coord_end < n_lines and lines[coord_end].strip() != "":
        coord_end += 1

    coordinates = lines[coord_start:coord_end]
    if not coordinates:
        raise ValueError("未找到坐标。")

    geometry = [charge_mult] + coordinates
    tail_start = coord_end + 1 if coord_end < n_lines else coord_end
    tail = lines[tail_start:]

    return title, geometry, tail


def add_block(output_parts, block):
    """Append a non-empty Gaussian job block followed by one blank separator line."""
    cleaned = block.rstrip("\n")
    if cleaned:
        output_parts.append(cleaned + "\n\n")


def process_one(
    path,
    args,
    header_template,
    scf_xqc_template,
    opt_link1_template,
    basis_template,
    stable_template,
    freq_template,
    sol_template,
    freq_only_template,
    sol_only_template,
    sol_from_freq_template,
):
    raw = (
        path.read_text(encoding="utf-8")
        .replace("\r\n", "\n")
        .replace("\r", "\n")
    )

    if args.check_link1 and (not args.force) and ("--Link1--" in raw):
        print(
            "[跳过] 已含 --Link1--：{} "
            "（使用 --no-check-link1 或 --force 可继续处理）".format(path)
        )
        return

    try:
        original_title, geometry, tail = split_sections(raw.split("\n"))
    except Exception as exc:
        print("[警告] 解析失败，跳过 {}：{}".format(path, exc))
        return

    if not geometry:
        print("[警告] 未找到电荷/多重度和坐标，跳过：{}".format(path))
        return

    charge_mult = geometry[0].strip()
    if not charge_mult:
        print("[警告] 电荷/多重度为空，跳过：{}".format(path))
        return

    base = path.stem
    use_genecp = args.basis_mode == "genecp"

    opt_based_job = args.jobs in {"all", "opt", "opt-freq", "opt-sol"}
    include_freq = args.jobs in {"all", "opt-freq", "freq", "freq-sol"}
    include_sol = args.jobs in {"all", "opt-sol", "sol", "freq-sol"}
    include_pre_stable = (
        opt_based_job
        and args.jobs != "opt"
        and (include_freq or include_sol)
        and (not args.skip_stable)
    )

    opt_keyword = (
        "opt"
        if args.opt_type == "minimum"
        else "opt=(ts,calcfc,noeigentest)"
    )

    dispersion_keyword = ""
    if args.dispersion.lower() not in {"none", "no", "off", ""}:
        dispersion_keyword = "EmpiricalDispersion={}".format(args.dispersion)

    basis_text = ""
    if use_genecp:
        basis_text = basis_template.format(
            nproc=args.nproc,
            mem=args.mem,
            solvent=args.solvent,
            functional=args.functional,
            dispersion_keyword=dispersion_keyword,
            charge_mult=charge_mult,
            base=base,
            genecp_feo_basis=args.genecp_feo_basis,
            genecp_chn_basis=args.genecp_chn_basis,
            genecp_feo_ecp=args.genecp_feo_ecp,
        ).strip()

    chk = "{}.chk".format(base)
    scf_chk = "{}_scf.chk".format(base)
    stable_chk = "{}_stable.chk".format(base)
    freq_chk = "{}_freq.chk".format(base)
    sol_chk = "{}_sol.chk".format(base)

    post_opt_chk = stable_chk if include_pre_stable else chk

    title_text = "\n".join(original_title).strip()
    if not title_text:
        title_text = "{} calculation".format(base)

    variables = {
        "base": base,
        "chk": chk,
        "scf_chk": scf_chk,
        "stable_chk": stable_chk,
        "freq_chk": freq_chk,
        "sol_chk": sol_chk,
        "post_opt_chk": post_opt_chk,
        "nproc": args.nproc,
        "mem": args.mem,
        "solvent": args.solvent,
        "sol_basis": args.sol_basis,
        "functional": args.functional,
        "dispersion_keyword": dispersion_keyword,
        "opt_keyword": opt_keyword,
        "opt_basis_keyword": (
            "genecp" if use_genecp else args.builtin_basis
        ),
        "stable_basis_keyword": (
            "genecp" if use_genecp else args.builtin_basis
        ),
        "freq_basis_keyword": (
            "genecp" if use_genecp else args.builtin_basis
        ),
        "charge_mult": charge_mult,
        "title": title_text,
        "opt_basis_block": basis_text,
        "stable_basis_block": basis_text,
        "freq_basis_block": basis_text,
        "sol_stable_keyword": (
            "" if args.no_sol_stable else "stable=opt"
        ),
    }

    new_header = header_template.format(**variables)
    new_scf_xqc = scf_xqc_template.format(**variables)
    new_opt_link1 = opt_link1_template.format(**variables)
    new_stable = stable_template.format(**variables)
    new_freq = freq_template.format(**variables)
    new_sol = sol_template.format(**variables)
    new_freq_only = freq_only_template.format(**variables)
    new_sol_only = sol_only_template.format(**variables)
    new_sol_from_freq = sol_from_freq_template.format(**variables)

    geometry_block = "\n".join(geometry).rstrip()
    tail_block = "\n".join(tail).rstrip() if tail else ""

    output_parts = []

    if opt_based_job:
        if args.pre_scf_xqc:
            # First job: converge the wavefunction with SCF=XQC.
            add_block(output_parts, new_scf_xqc)
            add_block(output_parts, geometry_block)

            if basis_text:
                add_block(output_parts, basis_text)

            # Second job: read the converged wavefunction and optimize.
            add_block(output_parts, new_opt_link1)
        else:
            # Default behavior: optimize directly.
            add_block(output_parts, new_header)
            add_block(output_parts, geometry_block)

            if basis_text:
                add_block(output_parts, basis_text)

        if include_pre_stable:
            add_block(output_parts, new_stable)

        if include_freq:
            add_block(output_parts, new_freq)

        if include_sol:
            add_block(output_parts, new_sol)

    elif args.jobs == "freq":
        # Pure frequency calculation: use the coordinates from the input directly.
        add_block(output_parts, new_freq_only)
        add_block(output_parts, geometry_block)
        if basis_text:
            add_block(output_parts, basis_text)

    elif args.jobs == "sol":
        # Pure SMD single point: use the coordinates from the input directly.
        add_block(output_parts, new_sol_only)
        add_block(output_parts, geometry_block)

    elif args.jobs == "freq-sol":
        # Pure Freq -> SMD: Freq starts from coordinates; SMD reads the Freq checkpoint.
        add_block(output_parts, new_freq_only)
        add_block(output_parts, geometry_block)
        if basis_text:
            add_block(output_parts, basis_text)
        add_block(output_parts, new_sol_from_freq)

    else:
        raise ValueError("未知 jobs 模式：{}".format(args.jobs))

    if args.keep_tail and tail_block:
        add_block(output_parts, tail_block)

    # One newline terminates the last content line; the next two create
    # two complete blank lines required at the end of the Gaussian input.
    new_text = "".join(output_parts).rstrip() + "\n\n\n"

    backup = path.with_suffix(path.suffix + ".bak")
    if not backup.exists():
        backup.write_text(raw, encoding="utf-8")

    path.write_text(new_text, encoding="utf-8")

    print(
        "[完成] {} | jobs={} | opt_type={} | pre_scf_xqc={} "
        "| pre_stable={} | sol_stable={} | check_link1={} | basis={} "
        "| functional={} | 备份={}".format(
            path,
            args.jobs,
            args.opt_type,
            args.pre_scf_xqc,
            include_pre_stable,
            not args.no_sol_stable,
            args.check_link1,
            args.basis_mode,
            args.functional,
            backup.name,
        )
    )



def prompt_text(question, default=None):
    """Read a text answer, returning the default when the user presses Enter."""
    suffix = " [{}]".format(default) if default is not None else ""
    answer = input("{}{}: ".format(question, suffix)).strip()
    if answer:
        return answer
    return default


def prompt_choice(question, choices, default_key):
    """
    Ask a numbered/keyed multiple-choice question.

    choices is an ordered list of (key, label, value).
    """
    print("\n{}".format(question))
    for key, label, _value in choices:
        default_mark = "（默认）" if key == default_key else ""
        print("  {}. {}{}".format(key, label, default_mark))

    valid = {key: value for key, _label, value in choices}

    while True:
        answer = input("请选择 [{}]: ".format(default_key)).strip()
        if not answer:
            answer = default_key
        if answer in valid:
            return valid[answer]
        print("输入无效，请选择：{}".format(
            ", ".join(key for key, _label, _value in choices)
        ))


def prompt_yes_no(question, default=False):
    """Ask a yes/no question and return a boolean."""
    default_text = "Y/n" if default else "y/N"
    while True:
        answer = input("{} [{}]: ".format(question, default_text)).strip().lower()
        if not answer:
            return default
        if answer in {"y", "yes", "是", "1"}:
            return True
        if answer in {"n", "no", "否", "0"}:
            return False
        print("请输入 y 或 n。")


def normalize_filename_keywords(values):
    """Return cleaned filename keywords, allowing repeated options or comma-separated input."""
    keywords = []
    for value in values or []:
        for item in value.split(","):
            cleaned = item.strip()
            if cleaned:
                keywords.append(cleaned)
    return keywords


def filename_matches_keywords(path, keywords, case_sensitive=False):
    """Return True when the filename contains at least one requested keyword."""
    if not keywords:
        return True

    filename = path.name
    if not case_sensitive:
        filename = filename.lower()
        keywords = [keyword.lower() for keyword in keywords]

    return any(keyword in filename for keyword in keywords)


def configure_interactively(args):
    """Populate argparse values through a Chinese question-and-answer flow."""
    print("\n=== Gaussian 输入批量生成向导 ===")
    print("核数固定为 32，内存固定为 60GB，不再询问。")

    args.target = prompt_text("待处理文件或目录", args.target or ".")

    interactive_target = Path(args.target).expanduser()
    if interactive_target.exists() and interactive_target.is_file():
        args.recursive = False
        print("检测到单个文件：仅处理 {}".format(interactive_target))
    else:
        args.recursive = prompt_yes_no("是否递归处理子文件夹", default=False)
        use_filename_filter = prompt_yes_no(
            "是否只处理文件名中含特定关键词的文件",
            default=False,
        )
        if use_filename_filter:
            keyword_text = prompt_text(
                "输入文件名关键词（多个关键词用英文逗号分隔；匹配任意一个即可）",
                "",
            )
            args.filename_keyword = normalize_filename_keywords([keyword_text])
            if not args.filename_keyword:
                print("未输入有效关键词，将处理所有匹配 --glob 的文件。")
            else:
                args.keyword_case_sensitive = prompt_yes_no(
                    "关键词匹配是否区分大小写",
                    default=False,
                )
        else:
            args.filename_keyword = []
            args.keyword_case_sensitive = False

    # 先选流程，这样纯 Freq / 纯 SMD 模式不会再询问无关的 Opt/TS 设置。
    args.jobs = prompt_choice(
        "选择计算流程",
        [
            ("1", "Opt → Stable=Opt → Freq → SMD", "all"),
            ("2", "仅 Opt", "opt"),
            ("3", "Opt → Stable=Opt → Freq", "opt-freq"),
            ("4", "Opt → Stable=Opt → SMD", "opt-sol"),
            ("5", "仅 Freq（不优化）", "freq"),
            ("6", "仅 SMD（不优化）", "sol"),
            ("7", "Freq → SMD（不优化）", "freq-sol"),
        ],
        default_key="1",
    )

    opt_based_job = args.jobs in {"all", "opt", "opt-freq", "opt-sol"}
    needs_gas_basis = args.jobs in {"all", "opt", "opt-freq", "opt-sol", "freq", "freq-sol"}
    needs_smd = args.jobs in {"all", "opt-sol", "sol", "freq-sol"}

    if opt_based_job:
        args.opt_type = prompt_choice(
            "选择优化类型",
            [
                ("1", "普通极小值优化（Opt）", "minimum"),
                ("2", "过渡态优化（Opt=TS）", "ts"),
            ],
            default_key="1",
        )
        args.pre_scf_xqc = prompt_yes_no(
            "是否先用 SCF=XQC 做第一段同类型优化",
            default=False,
        )
    else:
        args.opt_type = "minimum"
        args.pre_scf_xqc = False

    args.functional = prompt_text("泛函关键词", args.functional)

    args.dispersion = prompt_choice(
        "色散修正",
        [
            ("1", "不加色散", "none"),
            ("2", "D3(BJ)：EmpiricalDispersion=GD3BJ", "GD3BJ"),
            ("3", "D3：EmpiricalDispersion=GD3", "GD3"),
            ("4", "自定义色散关键词", "__custom__"),
        ],
        default_key="1",
    )
    if args.dispersion == "__custom__":
        args.dispersion = prompt_text("输入 EmpiricalDispersion 的值", "GD3BJ")

    if needs_gas_basis:
        args.basis_mode = prompt_choice(
            "Opt/Freq 基组模式",
            [
                (
                    "1",
                    "混合基组 GenECP（元素分组固定，基组名称自行设置）",
                    "genecp",
                ),
                ("2", "统一内置基组", "builtin"),
            ],
            default_key="1",
        )
        if args.basis_mode == "genecp":
            # 元素范围固定为 Fe/O 与 C/H/N，只允许修改对应的基组/ECP名称。
            args.genecp_feo_basis = prompt_text(
                "Fe/O 使用的轨道基组",
                args.genecp_feo_basis,
            )
            args.genecp_chn_basis = prompt_text(
                "C/H/N 使用的轨道基组",
                args.genecp_chn_basis,
            )
            args.genecp_feo_ecp = prompt_text(
                "Fe/O 使用的 ECP 名称",
                args.genecp_feo_ecp,
            )
        else:
            args.builtin_basis = prompt_text(
                "Opt/Freq 使用的统一基组",
                args.builtin_basis,
            )

    args.skip_stable = False
    if args.jobs in {"all", "opt-freq", "opt-sol"}:
        args.skip_stable = not prompt_yes_no(
            "Opt 后、Freq/SMD 前是否执行气相 Stable=Opt",
            default=True,
        )

    if needs_smd:
        args.solvent = prompt_text("SMD 溶剂", args.solvent)
        args.sol_basis = prompt_text("SMD 单点基组", args.sol_basis)
        args.no_sol_stable = not prompt_yes_no(
            "SMD 单点是否加入 Stable=Opt",
            default=True,
        )
    else:
        args.no_sol_stable = False

    args.force = prompt_yes_no(
        "遇到已经含有 --Link1-- 的文件时是否强制重写",
        default=False,
    )
    if args.force:
        args.check_link1 = False

    print("\n=== 本次设置 ===")
    target_display = Path(args.target).expanduser()
    print("目标：{}".format(target_display))
    if target_display.exists() and target_display.is_file():
        print("处理模式：单个文件")
    else:
        print("处理模式：目录")
        print("递归：{}".format("是" if args.recursive else "否"))
        normalized_keywords = normalize_filename_keywords(args.filename_keyword)
        if normalized_keywords:
            print("文件名关键词：{}".format("，".join(normalized_keywords)))
            print("关键词区分大小写：{}".format(
                "是" if args.keyword_case_sensitive else "否"
            ))
        else:
            print("文件名关键词：不筛选")

    print("流程：{}".format(args.jobs))
    if opt_based_job:
        print("优化：{}".format(
            "TS" if args.opt_type == "ts" else "普通极小值"
        ))
        print("第一段 SCF=XQC：{}".format(
            "是" if args.pre_scf_xqc else "否"
        ))
    print("泛函：{}".format(args.functional))
    print("色散：{}".format(args.dispersion))

    if needs_gas_basis:
        print("Opt/Freq 基组模式：{}".format(args.basis_mode))
        if args.basis_mode == "genecp":
            print("GenECP Fe/O 轨道基组：{}".format(args.genecp_feo_basis))
            print("GenECP C/H/N 轨道基组：{}".format(args.genecp_chn_basis))
            print("GenECP Fe/O ECP：{}".format(args.genecp_feo_ecp))
        else:
            print("Opt/Freq 统一基组：{}".format(args.builtin_basis))

    print("气相 Stable=Opt：{}".format(
        "是"
        if args.jobs in {"all", "opt-freq", "opt-sol"} and not args.skip_stable
        else "否"
    ))
    if needs_smd:
        print("SMD：{} / {}".format(args.solvent, args.sol_basis))
        print("SMD Stable=Opt：{}".format(
            "否" if args.no_sol_stable else "是"
        ))
    print("核数/内存：32 / 60GB")

    if not prompt_yes_no("\n确认生成输入文件", default=True):
        print("已取消。")
        sys.exit(0)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "批量重写 Gaussian .gjf：支持 Opt/Stable/Freq/SMD 组合，"
            "也支持纯 Freq、纯 SMD、Freq→SMD（不优化）"
        )
    )

    parser.add_argument(
        "--interactive",
        action="store_true",
        help="进入问答式配置；不带任何参数运行时也会自动进入",
    )

    parser.add_argument(
        "target",
        nargs="?",
        default=".",
        help=(
            "待处理的单个 .gjf 文件或目录；默认当前目录。"
            "例如：python process_interactive.py 2int2.gjf，"
            "或 python process_interactive.py symmL/"
        ),
    )

    parser.add_argument("-r", "--recursive", action="store_true")
    parser.add_argument(
        "--filename-keyword",
        action="append",
        default=[],
        help=(
            "目录模式下仅处理文件名包含该关键词的文件；可重复使用，"
            "也可用英文逗号分隔多个关键词；多个关键词之间为任意匹配"
        ),
    )
    parser.add_argument(
        "--keyword-case-sensitive",
        action="store_true",
        help="文件名关键词匹配区分大小写；默认不区分",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="强制处理已经含有 --Link1-- 的文件；保留用于兼容旧命令",
    )

    link1_check_group = parser.add_mutually_exclusive_group()
    link1_check_group.add_argument(
        "--check-link1",
        dest="check_link1",
        action="store_true",
        default=True,
        help="处理前检查文件是否已有 --Link1--；这是默认行为",
    )
    link1_check_group.add_argument(
        "--no-check-link1",
        dest="check_link1",
        action="store_false",
        help="关闭已有 --Link1-- 检查，仍然继续重写文件",
    )

    parser.add_argument(
        "--keep-tail",
        action="store_true",
        help="保留原坐标后的尾部内容；默认丢弃，防止旧基组重复写入",
    )
    parser.add_argument(
        "--pre-scf-xqc",
        action="store_true",
        help=(
            "先用 SCF=XQC 做一次同类型优化，再通过 Link1 继续同类型优化；"
            "默认关闭"
        ),
    )
    parser.add_argument(
        "--skip-stable",
        action="store_true",
        help="跳过 Opt 后、Freq/SMD 前的气相 Stable=Opt",
    )
    parser.add_argument(
        "--no-sol-stable",
        action="store_true",
        help="溶剂单点中不加入 Stable=Opt",
    )

    parser.add_argument("--nproc", default="32")
    parser.add_argument("--mem", default="60GB")
    parser.add_argument("--solvent", default="methanol")
    parser.add_argument("--sol-basis", default="def2tzvp")

    parser.add_argument(
        "--functional",
        default="UPW6B95",
        help="泛函关键词；默认 UPW6B95",
    )
    parser.add_argument(
        "--dispersion",
        default="GD3BJ",
        help="EmpiricalDispersion 参数；默认 GD3BJ；写 none 可关闭",
    )

    parser.add_argument(
        "--jobs",
        choices=["all", "opt", "opt-freq", "opt-sol", "freq", "sol", "freq-sol"],
        default="all",
        help=(
            "all=Opt+Stable+Freq+SMD；"
            "opt=仅 Opt；"
            "opt-freq=Opt+Stable+Freq；"
            "opt-sol=Opt+Stable+SMD；"
            "freq=仅 Freq；"
            "sol=仅 SMD；"
            "freq-sol=Freq+SMD（不优化）"
        ),
    )

    parser.add_argument(
        "--opt-type",
        choices=["minimum", "ts"],
        default="minimum",
        help="minimum=普通极小值；ts=过渡态",
    )

    parser.add_argument(
        "--basis-mode",
        choices=["genecp", "builtin"],
        default="genecp",
        help="genecp=写混合基组块；builtin=使用统一内置基组",
    )
    parser.add_argument(
        "--builtin-basis",
        default="def2svp",
        help="basis-mode=builtin 时，Opt/Stable/Freq 使用的基组",
    )
    parser.add_argument(
        "--genecp-feo-basis",
        default="def2tzvp",
        help="basis-mode=genecp 时 Fe/O 使用的轨道基组；元素范围固定为 Fe O",
    )
    parser.add_argument(
        "--genecp-chn-basis",
        default="def2svp",
        help="basis-mode=genecp 时 C/H/N 使用的轨道基组；元素范围固定为 C H N",
    )
    parser.add_argument(
        "--genecp-feo-ecp",
        default="def2tzvp",
        help="basis-mode=genecp 时 Fe/O 使用的 ECP 名称；元素范围固定为 Fe O",
    )

    parser.add_argument("--header", default=None)
    parser.add_argument("--scf-xqc", default=None, help="预 SCF=XQC 模板")
    parser.add_argument(
        "--opt-link1",
        default=None,
        help="SCF=XQC 后的优化 Link1 模板",
    )
    parser.add_argument("--basis", default=None)
    parser.add_argument("--stable", default=None)
    parser.add_argument("--link1a", default=None, help="频率 Link1 模板")
    parser.add_argument("--link1b", default=None, help="溶剂单点 Link1 模板")
    parser.add_argument("--glob", default="*.gjf")

    no_cli_options = len(sys.argv) == 1
    args = parser.parse_args()

    if args.interactive or no_cli_options:
        configure_interactively(args)

    try:
        header_template = read_template(args.header, DEFAULT_HEADER)
        scf_xqc_template = read_template(args.scf_xqc, DEFAULT_SCF_XQC)
        opt_link1_template = read_template(
            args.opt_link1,
            DEFAULT_OPT_LINK1,
        )
        basis_template = read_template(args.basis, DEFAULT_BASIS)
        stable_template = read_template(args.stable, DEFAULT_STABLE)
        freq_template = read_template(args.link1a, DEFAULT_FREQ)
        sol_template = read_template(args.link1b, DEFAULT_SOL)
        freq_only_template = DEFAULT_FREQ_ONLY
        sol_only_template = DEFAULT_SOL_ONLY
        sol_from_freq_template = DEFAULT_SOL_FROM_FREQ
    except Exception as exc:
        print("[错误] 读取模板失败：{}".format(exc))
        sys.exit(1)

    target = Path(args.target).expanduser().resolve()

    if not target.exists():
        print("[错误] 目标路径不存在：{}".format(target))
        sys.exit(1)

    if target.is_file():
        # Supplying a file path always means: process exactly this one file.
        # The --glob and --recursive options are intentionally ignored.
        files = [target]
        print("[单文件模式] {}".format(target))
    elif target.is_dir():
        candidate_files = (
            list(target.rglob(args.glob))
            if args.recursive
            else list(target.glob(args.glob))
        )
        candidate_files = sorted(
            path for path in candidate_files if path.is_file()
        )

        if not candidate_files:
            print("在 {} 中未找到匹配 {} 的文件".format(target, args.glob))
            sys.exit(0)

        filename_keywords = normalize_filename_keywords(args.filename_keyword)
        files = [
            path
            for path in candidate_files
            if filename_matches_keywords(
                path,
                filename_keywords,
                case_sensitive=args.keyword_case_sensitive,
            )
        ]

        if filename_keywords:
            print(
                "[文件名筛选] 关键词={} | 区分大小写={} | 候选={} | 命中={}".format(
                    ", ".join(filename_keywords),
                    args.keyword_case_sensitive,
                    len(candidate_files),
                    len(files),
                )
            )

        if not files:
            print(
                "在 {} 中找到 {} 个匹配 {} 的文件，但没有文件名包含关键词：{}".format(
                    target,
                    len(candidate_files),
                    args.glob,
                    ", ".join(filename_keywords),
                )
            )
            sys.exit(0)
    else:
        print("[错误] 目标路径既不是普通文件也不是目录：{}".format(target))
        sys.exit(1)

    for path in files:
        try:
            process_one(
                path,
                args,
                header_template,
                scf_xqc_template,
                opt_link1_template,
                basis_template,
                stable_template,
                freq_template,
                sol_template,
                freq_only_template,
                sol_only_template,
                sol_from_freq_template,
            )
        except Exception as exc:
            print("[错误] {}：{}".format(path, exc))


if __name__ == "__main__":
    main()









