#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
通用交互式 Gaussian/GJF 批量关键词替换工具。

特点：
1. 不预设 B3LYP 或任何目标泛函/关键词。
2. 用户自行输入：查找内容、替换内容、匹配范围和匹配方式。
3. 可一次输入多个替换值，为每个替换值生成一套新文件。
4. 可替换文件名、Gaussian route section、整个文件内容，或组合替换。
5. 默认保留原文件；也支持原地修改并自动创建 .bak 备份。
6. 支持单文件、文件夹、递归子文件夹和自定义 glob。
7. 可在替换的同时，通过命令行或交互模式向 route section 或文件末尾追加内容。

直接运行：
    python interactive_batch_replace_gjf.py

也可指定目标后进入交互：
    python interactive_batch_replace_gjf.py ./gjf_folder
"""

import argparse
import re
import shutil
import sys
from pathlib import Path
from typing import Callable, Iterable, List, Optional, Sequence, Tuple, Pattern


def prompt_text(question: str, default: Optional[str] = None, allow_empty: bool = False) -> str:
    suffix = f" [{default}]" if default is not None else ""
    while True:
        answer = input(f"{question}{suffix}: ").strip()
        if answer:
            return answer
        if default is not None:
            return default
        if allow_empty:
            return ""
        print("输入不能为空。")


def prompt_yes_no(question: str, default: bool = False) -> bool:
    marker = "Y/n" if default else "y/N"
    while True:
        answer = input(f"{question} [{marker}]: ").strip().lower()
        if not answer:
            return default
        if answer in {"y", "yes", "1", "是"}:
            return True
        if answer in {"n", "no", "0", "否"}:
            return False
        print("请输入 y 或 n。")


def prompt_choice(
    question: str,
    choices: Sequence[Tuple[str, str, str]],
    default_key: str,
) -> str:
    print(f"\n{question}")
    values = {}
    for key, label, value in choices:
        default_mark = "（默认）" if key == default_key else ""
        print(f"  {key}. {label}{default_mark}")
        values[key] = value

    while True:
        answer = input(f"请选择 [{default_key}]: ").strip() or default_key
        if answer in values:
            return values[answer]
        print("输入无效，请重新选择。")


def split_replacements(raw: str) -> List[str]:
    """
    支持用英文逗号、中文逗号或分号分隔多个替换值。
    若替换文本本身含逗号，可在交互中选择逐行输入模式。
    """
    values = [item.strip() for item in re.split(r"[,，;；]", raw)]
    return [item for item in values if item]


def read_replacements() -> List[str]:
    mode = prompt_choice(
        "如何输入替换值",
        [
            ("1", "输入一个替换值", "single"),
            ("2", "一行输入多个值，用逗号或分号分隔", "separated"),
            ("3", "逐行输入多个值，直接回车结束", "lines"),
        ],
        default_key="1",
    )

    if mode == "single":
        return [prompt_text("替换成什么")]

    if mode == "separated":
        while True:
            values = split_replacements(prompt_text("输入替换值列表"))
            if values:
                return values
            print("至少需要一个替换值。")

    values: List[str] = []
    print("逐行输入替换值；输入空行结束：")
    while True:
        value = input(f"  替换值 {len(values) + 1}: ").strip()
        if not value:
            if values:
                return values
            print("至少需要一个替换值。")
            continue
        values.append(value)


def build_pattern(old: str, match_mode: str, case_sensitive: bool) -> Pattern[str]:
    flags = 0 if case_sensitive else re.IGNORECASE

    if match_mode == "literal":
        expression = re.escape(old)
    elif match_mode == "keyword":
        # 只把字母和数字视为关键词内部字符；下划线和连字符可作为分隔符。
        # 因此 test_b3lyp.gjf、B3LYP/、M06-2X/ 等都可正常匹配。
        expression = rf"(?<![A-Za-z0-9]){re.escape(old)}(?![A-Za-z0-9])"
    elif match_mode == "regex":
        expression = old
    else:
        raise ValueError(f"未知匹配模式：{match_mode}")

    try:
        return re.compile(expression, flags)
    except re.error as exc:
        raise ValueError(f"正则表达式无效：{exc}") from exc


def replace_text(text: str, pattern: Pattern[str], replacement: str) -> Tuple[str, int]:
    """正则模式下允许 \1 等反向引用；普通模式下 replacement 按字面替换。"""
    if pattern.pattern and "\\" in replacement:
        try:
            return pattern.subn(replacement, text)
        except re.error:
            # 对普通替换词中的反斜杠提供安全回退。
            return pattern.subn(lambda _match: replacement, text)
    return pattern.subn(lambda _match: replacement, text)


def find_route_range(lines: Sequence[str]) -> Optional[Tuple[int, int]]:
    """返回 Gaussian route section 的 [start, end) 行号范围。"""
    start = None
    for index, line in enumerate(lines):
        if line.lstrip().startswith("#"):
            start = index
            break

    if start is None:
        return None

    end = len(lines)
    for index in range(start + 1, len(lines)):
        if lines[index].strip() == "":
            end = index
            break

    return start, end


def replace_in_route(
    text: str,
    pattern: Pattern[str],
    replacement: str,
) -> Tuple[str, int]:
    lines = text.splitlines(keepends=True)
    route_range = find_route_range(lines)
    if route_range is None:
        return text, 0

    start, end = route_range
    total = 0
    for index in range(start, end):
        lines[index], count = replace_text(lines[index], pattern, replacement)
        total += count

    return "".join(lines), total


def detect_newline(text: str) -> str:
    """尽量沿用原文件换行符。"""
    if "\r\n" in text:
        return "\r\n"
    if "\r" in text:
        return "\r"
    return "\n"


def resolve_added_content(value: str) -> str:
    """
    解析额外内容。

    普通字符串直接使用；若以 @ 开头，则读取其后的文本文件。
    例如：--append-end @pcm_parameters.txt
    """
    if not value.startswith("@"):
        return value

    path = Path(value[1:]).expanduser()
    if not path.exists() or not path.is_file():
        raise ValueError(f"额外内容文件不存在：{path}")

    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"额外内容文件不是 UTF-8 编码：{path}") from exc
    except OSError as exc:
        raise ValueError(f"无法读取额外内容文件：{path}：{exc}") from exc


def add_to_route(text: str, additions: Sequence[str]) -> Tuple[str, int]:
    """
    将额外关键词/文本作为 route section 的续行，插入到 route 后的空行之前。

    例如：
        #p PBE0/def2svp
        scf=xqc int=ultrafine

    Gaussian 允许 route section 跨多行。
    """
    blocks = [block.strip("\r\n") for block in additions if block.strip()]
    if not blocks:
        return text, 0

    lines = text.splitlines(keepends=True)
    route_range = find_route_range(lines)
    if route_range is None:
        return text, 0

    _start, end = route_range
    newline = detect_newline(text)
    inserted_lines: List[str] = []

    for block in blocks:
        normalized = block.replace("\r\n", "\n").replace("\r", "\n")
        for line in normalized.split("\n"):
            if line.strip():
                inserted_lines.append(line.rstrip() + newline)

    if not inserted_lines:
        return text, 0

    lines[end:end] = inserted_lines
    return "".join(lines), len(blocks)


def append_to_file_end(text: str, additions: Sequence[str]) -> Tuple[str, int]:
    """把额外内容作为独立输入块追加到文件末尾，并以空行结束。"""
    blocks = [block.strip("\r\n") for block in additions if block.strip()]
    if not blocks:
        return text, 0

    newline = detect_newline(text)
    base = text.rstrip("\r\n")
    normalized_blocks = []

    for block in blocks:
        normalized = block.replace("\r\n", "\n").replace("\r", "\n")
        normalized_blocks.append(normalized.replace("\n", newline).rstrip("\r\n"))

    result = base + newline * 2 + (newline * 2).join(normalized_blocks) + newline * 2
    return result, len(normalized_blocks)


def safe_filename_replacement(replacement: str) -> str:
    """替换 Windows/Linux 文件名中不安全的字符。"""
    safe = replacement.strip()
    safe = re.sub(r"[\\/:*?\"<>|]", "_", safe)
    safe = re.sub(r"\s+", "_", safe)
    return safe or "replacement"


def replace_filename(
    path: Path,
    pattern: Pattern[str],
    replacement: str,
) -> Tuple[str, int]:
    safe_replacement = safe_filename_replacement(replacement)
    return replace_text(path.name, pattern, safe_replacement)


def collect_files(
    target: Path,
    recursive: bool,
    glob_pattern: str,
    filename_filter: Optional[Pattern[str]],
) -> List[Path]:
    if target.is_file():
        return [target]

    iterator: Iterable[Path]
    iterator = target.rglob(glob_pattern) if recursive else target.glob(glob_pattern)

    files = [path for path in iterator if path.is_file()]
    if filename_filter is not None:
        files = [path for path in files if filename_filter.search(path.name)]

    return sorted(files)


def make_copy_destination(
    source: Path,
    output_name: str,
    replacement: str,
    name_was_replaced: bool,
) -> Path:
    if name_was_replaced:
        return source.with_name(output_name)

    suffix = safe_filename_replacement(replacement)
    return source.with_name(f"{source.stem}_{suffix}{source.suffix}")


def configure_interactively(args: argparse.Namespace) -> argparse.Namespace:
    print("\n=== 通用 GJF 批量关键词替换向导 ===")

    args.target = prompt_text("待处理的单个文件或文件夹", args.target or ".")
    target = Path(args.target).expanduser()

    if target.exists() and target.is_file():
        args.recursive = False
        print(f"检测到单个文件：{target}")
    else:
        args.recursive = prompt_yes_no("是否递归处理子文件夹", default=False)

    args.glob = prompt_text("文件匹配规则", args.glob or "*.gjf")

    args.old = prompt_text("要查找/替换的旧关键词或文本")
    args.replacements = read_replacements()

    args.scope = prompt_choice(
        "替换范围",
        [
            ("1", "文件名 + Gaussian route section（# 开头部分）", "name_route"),
            ("2", "仅 Gaussian route section", "route"),
            ("3", "文件名 + 整个文件内容", "name_all"),
            ("4", "仅整个文件内容", "all"),
            ("5", "仅文件名", "name"),
        ],
        default_key="1",
    )

    args.match_mode = prompt_choice(
        "匹配方式",
        [
            ("1", "普通文本/子字符串匹配", "literal"),
            ("2", "独立关键词匹配，避免改到更长单词内部", "keyword"),
            ("3", "正则表达式（高级）", "regex"),
        ],
        default_key="2",
    )

    args.case_sensitive = prompt_yes_no("是否区分大小写", default=False)

    args.filter_names = prompt_yes_no(
        "是否只处理文件名中含旧关键词的文件",
        default=True,
    )

    args.in_place = prompt_yes_no(
        "是否直接修改原文件（否则生成新副本）",
        default=False,
    )

    if args.in_place:
        if len(args.replacements) > 1:
            print("[提示] 原地修改只能使用一个替换值，将使用第一个。")
            args.replacements = args.replacements[:1]
        args.backup = prompt_yes_no("原地修改前是否创建 .bak 备份", default=True)
        args.overwrite = True
    else:
        args.backup = False
        args.overwrite = prompt_yes_no("目标副本已存在时是否覆盖", default=False)

    args.add_route = []
    args.append_end = []
    if prompt_yes_no("是否还要额外增加内容", default=False):
        print("提示：输入 @文件名 可从文本文件读取多行内容。")
        while True:
            add_position = prompt_choice(
                "增加到什么位置",
                [
                    ("1", "追加到 Gaussian route section（#p 区域）", "route"),
                    ("2", "作为独立参数块追加到文件末尾", "end"),
                    ("3", "完成，不再增加", "done"),
                ],
                default_key="3",
            )
            if add_position == "done":
                break

            value = prompt_text("要增加的内容，或输入 @文本文件路径")
            if add_position == "route":
                args.add_route.append(value)
            else:
                args.append_end.append(value)

    print("\n=== 本次设置 ===")
    print(f"目标：{Path(args.target).expanduser()}")
    print(f"匹配文件：{args.glob}")
    print(f"递归：{'是' if args.recursive else '否'}")
    print(f"旧文本：{args.old}")
    print(f"替换值：{', '.join(args.replacements)}")
    print(f"替换范围：{args.scope}")
    print(f"匹配模式：{args.match_mode}")
    print(f"区分大小写：{'是' if args.case_sensitive else '否'}")
    print(f"文件名预筛选：{'是' if args.filter_names else '否'}")
    print(f"输出方式：{'原地修改' if args.in_place else '生成副本'}")
    if args.add_route:
        print(f"route 追加内容：{len(args.add_route)} 项")
    if args.append_end:
        print(f"文件末尾追加内容：{len(args.append_end)} 项")

    if not prompt_yes_no("\n确认开始处理", default=True):
        print("已取消。")
        sys.exit(0)

    return args


def process_one(
    source: Path,
    pattern: Pattern[str],
    replacements: Sequence[str],
    scope: str,
    in_place: bool,
    backup: bool,
    overwrite: bool,
    encoding: str,
    add_route_blocks: Sequence[str],
    append_end_blocks: Sequence[str],
) -> Tuple[int, int]:
    try:
        original = source.read_text(encoding=encoding)
    except UnicodeDecodeError:
        print(f"[错误] 无法使用 {encoding} 解码：{source}")
        return 0, 1
    except OSError as exc:
        print(f"[错误] 读取失败：{source}：{exc}")
        return 0, 1

    generated = 0
    skipped = 0

    for replacement in replacements:
        new_text = original
        content_count = 0
        filename_count = 0
        output_name = source.name

        if scope in {"route", "name_route"}:
            new_text, content_count = replace_in_route(new_text, pattern, replacement)
        elif scope in {"all", "name_all"}:
            new_text, content_count = replace_text(new_text, pattern, replacement)

        if scope in {"name", "name_route", "name_all"}:
            output_name, filename_count = replace_filename(source, pattern, replacement)

        addition_count = 0
        if add_route_blocks:
            new_text, count = add_to_route(new_text, add_route_blocks)
            addition_count += count
        if append_end_blocks:
            new_text, count = append_to_file_end(new_text, append_end_blocks)
            addition_count += count

        if content_count == 0 and filename_count == 0 and addition_count == 0:
            print(f"[跳过] 未找到匹配，且没有成功增加内容：{source}")
            skipped += 1
            continue

        if in_place:
            destination = source
            if backup:
                backup_path = source.with_suffix(source.suffix + ".bak")
                if not backup_path.exists():
                    try:
                        shutil.copy2(source, backup_path)
                    except OSError as exc:
                        print(f"[错误] 备份失败：{source}：{exc}")
                        skipped += 1
                        continue

            try:
                destination.write_text(new_text, encoding=encoding)
            except OSError as exc:
                print(f"[错误] 写入失败：{destination}：{exc}")
                skipped += 1
                continue

            # 原地模式下如果选择替换文件名，则写完内容后再重命名。
            if filename_count > 0 and output_name != source.name:
                renamed = source.with_name(output_name)
                if renamed.exists() and renamed != source:
                    print(f"[错误] 重命名目标已存在：{renamed}")
                    skipped += 1
                    continue
                try:
                    source.rename(renamed)
                    destination = renamed
                except OSError as exc:
                    print(f"[错误] 重命名失败：{source}：{exc}")
                    skipped += 1
                    continue

            print(
                f"[修改] {source.name} -> {destination.name} | "
                f"内容替换={content_count} | 文件名替换={filename_count} | "
                f"新增内容块={addition_count}"
            )
            generated += 1
            continue

        destination = make_copy_destination(
            source=source,
            output_name=output_name,
            replacement=replacement,
            name_was_replaced=(filename_count > 0),
        )

        if destination.resolve() == source.resolve():
            destination = source.with_name(
                f"{source.stem}_{safe_filename_replacement(replacement)}{source.suffix}"
            )

        if destination.exists() and not overwrite:
            print(f"[跳过] 目标文件已存在：{destination}")
            skipped += 1
            continue

        try:
            destination.write_text(new_text, encoding=encoding)
        except OSError as exc:
            print(f"[错误] 写入失败：{destination}：{exc}")
            skipped += 1
            continue

        print(
            f"[生成] {source.name} -> {destination.name} | "
            f"替换值={replacement} | 内容替换={content_count} | "
            f"文件名替换={filename_count} | 新增内容块={addition_count}"
        )
        generated += 1

    return generated, skipped


def main() -> None:
    parser = argparse.ArgumentParser(
        description="通用交互式 Gaussian/GJF 批量关键词替换工具"
    )
    parser.add_argument("target", nargs="?", default=".")
    parser.add_argument("--interactive", action="store_true")
    parser.add_argument("--glob", default="*.gjf")
    parser.add_argument("-r", "--recursive", action="store_true")
    parser.add_argument("--encoding", default="utf-8")

    # 以下参数也保留命令行用法；不完整时自动进入交互模式。
    parser.add_argument("--old", default=None, help="要查找的旧文本")
    parser.add_argument(
        "--replace",
        dest="replacements",
        action="append",
        help="替换值；可重复使用多次",
    )
    parser.add_argument(
        "--scope",
        choices=["name_route", "route", "name_all", "all", "name"],
        default=None,
    )
    parser.add_argument(
        "--match-mode",
        choices=["literal", "keyword", "regex"],
        default=None,
    )
    parser.add_argument("--case-sensitive", action="store_true")
    parser.add_argument("--no-name-filter", dest="filter_names", action="store_false")
    parser.set_defaults(filter_names=True)
    parser.add_argument("--in-place", action="store_true")
    parser.add_argument("--no-backup", dest="backup", action="store_false")
    parser.set_defaults(backup=True)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--add-route",
        action="append",
        default=[],
        metavar="TEXT",
        help=(
            "向 Gaussian route section 追加内容；可重复使用。"
            "TEXT 以 @ 开头时，从 UTF-8 文本文件读取"
        ),
    )
    parser.add_argument(
        "--append-end",
        action="append",
        default=[],
        metavar="TEXT",
        help=(
            "把内容作为独立参数块追加到文件末尾；可重复使用。"
            "TEXT 以 @ 开头时，从 UTF-8 文本文件读取"
        ),
    )

    args = parser.parse_args()

    incomplete_cli = not (
        args.old is not None
        and args.replacements
        and args.scope is not None
        and args.match_mode is not None
    )

    if args.interactive or len(sys.argv) == 1 or incomplete_cli:
        configure_interactively(args)

    target = Path(args.target).expanduser().resolve()
    if not target.exists():
        print(f"[错误] 目标路径不存在：{target}")
        sys.exit(1)

    if args.in_place and len(args.replacements) > 1:
        print("[错误] 原地修改模式只能使用一个替换值。")
        sys.exit(1)

    try:
        add_route_blocks = [resolve_added_content(value) for value in args.add_route]
        append_end_blocks = [resolve_added_content(value) for value in args.append_end]
    except ValueError as exc:
        print(f"[错误] {exc}")
        sys.exit(1)

    try:
        pattern = build_pattern(args.old, args.match_mode, args.case_sensitive)
    except ValueError as exc:
        print(f"[错误] {exc}")
        sys.exit(1)

    filename_filter = pattern if args.filter_names and target.is_dir() else None
    files = collect_files(
        target=target,
        recursive=args.recursive,
        glob_pattern=args.glob,
        filename_filter=filename_filter,
    )

    if not files:
        extra = "且文件名匹配旧关键词" if filename_filter is not None else ""
        print(f"[结束] 未找到匹配 {args.glob} {extra}的文件：{target}")
        return

    print(f"\n[找到] {len(files)} 个待处理文件")

    total_generated = 0
    total_skipped = 0

    for source in files:
        generated, skipped = process_one(
            source=source,
            pattern=pattern,
            replacements=args.replacements,
            scope=args.scope,
            in_place=args.in_place,
            backup=args.backup,
            overwrite=args.overwrite,
            encoding=args.encoding,
            add_route_blocks=add_route_blocks,
            append_end_blocks=append_end_blocks,
        )
        total_generated += generated
        total_skipped += skipped

    print(
        f"\n[完成] 输入文件={len(files)} | "
        f"成功输出/修改={total_generated} | 跳过/失败={total_skipped}"
    )


if __name__ == "__main__":
    main()







