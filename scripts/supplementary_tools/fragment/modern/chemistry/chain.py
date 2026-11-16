# Shared molecular chain traversal helpers.
from .constants import SINGLE_CONNECT_ATOM_OR_GROUP as SCAOG


def search_chain(terminal_atom_arr, connect_dict, additional_conditions=None):
    if len(connect_dict) == 1:
        return [list(connect_dict.keys())]
    length_between_terminal = length_two_atoms(terminal_atom_arr, connect_dict)

    main_chain = get_lengthest_chain(
        length_between_terminal, additional_conditions=additional_conditions
    )
    all_chain = [length_between_terminal[main_chain]]
    branch_links(
        main_chain,
        length_between_terminal,
        connect_dict,
        terminal_atom_arr,
        additional_conditions=additional_conditions,
        all_chain=all_chain,
    )  # noqa
    return all_chain


def subchain_find(main_chain, length_between_terminal, j, k):
    start_atom = int(main_chain.split("-")[0])
    sub_chain_of_this_atom = {}
    sub_side_chain_name = ""
    for item in length_between_terminal.keys():
        if item != main_chain:
            string = item.split("-")
            if int(string[0]) == start_atom or int(string[1]) == start_atom:
                if (
                    j in length_between_terminal[item]
                    and k in length_between_terminal[item]
                ):
                    a = length_between_terminal[item].index(j)
                    b = length_between_terminal[item].index(k)
                    if a < b:
                        side_arr = length_between_terminal[item][
                            a : len(length_between_terminal[item])
                        ]
                    elif a > b:
                        side_arr = length_between_terminal[item][0 : a + 1]
                        side_arr.reverse()
                    sub_side_chain_name = str(j) + "-" + str(side_arr[-1])
                    sub_chain_of_this_atom[sub_side_chain_name] = side_arr
    return (sub_chain_of_this_atom, sub_side_chain_name)


def get_lengthest_chain(length_between_terminal, additional_conditions=None):
    if additional_conditions is not None:
        nn = 0
        tmp = {}
        for item, atoms in length_between_terminal.items():
            n = 0
            for ac in atoms:
                n += additional_conditions[ac]
            tmp[item] = n
            if n > nn:
                nn = n
        longest_chain_arr = [item for item in tmp.keys() if tmp[item] == nn]
    else:
        longest_chain_arr = list(length_between_terminal.keys())
    if len(longest_chain_arr) == 1:
        return longest_chain_arr[0]
    tmp = {}
    nn = 0
    for item in longest_chain_arr:
        tmp[item] = len(length_between_terminal[item])
        if len(length_between_terminal[item]) > nn:
            nn = len(length_between_terminal[item])
    longest_chain_arr = [item for item in tmp.keys() if tmp[item] == nn]
    return longest_chain_arr[0]


def branch_links(
    main_chain,
    length_between_terminal,
    connect_dict,
    terminal_arr,
    additional_conditions=None,
    *,
    all_chain,
):
    main_chain_n = len(length_between_terminal[main_chain])
    for i in range(1, main_chain_n - 1):
        j = length_between_terminal[main_chain][i]
        if (
            len(
                set(connect_dict[j]).difference(
                    set(length_between_terminal[main_chain])
                )
            )
            > 0
        ):
            for k in connect_dict[j]:
                if (
                    k not in length_between_terminal[main_chain]
                    and k not in terminal_arr
                ):
                    sub_side_chain, sub_side_chain_name = subchain_find(
                        main_chain, length_between_terminal, j, k
                    )  # noqa
                    if len(sub_side_chain) == 1:
                        all_chain.append(sub_side_chain[sub_side_chain_name])

                    elif len(sub_side_chain) == 0:
                        pass
                    else:
                        main_chain_this_side_chain = get_lengthest_chain(
                            sub_side_chain, additional_conditions=additional_conditions
                        )  # noqa
                        all_chain.append(sub_side_chain[main_chain_this_side_chain])
                        branch_links(
                            main_chain_this_side_chain,
                            sub_side_chain,
                            connect_dict,
                            terminal_arr,
                            additional_conditions=additional_conditions,
                            all_chain=all_chain,
                        )


def length_two_atoms(terminal_atom_arr, reduce_connect):
    terminal_n = len(terminal_atom_arr)
    length_between_terminal = {}
    for i in range(0, terminal_n):
        j = terminal_atom_arr[i]
        break_flag = 0
        chain_atom_arr = [terminal_atom_arr[i]]
        del_atom = []
        while 1:
            now_center_atom = j
            for k in reduce_connect[now_center_atom]:
                if k not in chain_atom_arr and k not in del_atom:
                    if k in terminal_atom_arr:
                        if k > terminal_atom_arr[i]:
                            chain_name = str(terminal_atom_arr[i]) + "-" + str(k)
                        else:
                            chain_name = str(k) + "-" + str(terminal_atom_arr[i])
                        chain_atom_arr.append(k)
                        if chain_name not in length_between_terminal.keys():
                            length_between_terminal[chain_name] = chain_atom_arr[:]
                        break_flag += 1
                        break_flag_inner = 0
                        while 1:
                            if chain_atom_arr[-1] == terminal_atom_arr[i]:
                                j = terminal_atom_arr[i]
                                break
                            else:
                                back_atom = chain_atom_arr[-1]
                                for kk in reduce_connect[back_atom]:
                                    if kk not in chain_atom_arr and kk not in del_atom:
                                        if kk in terminal_atom_arr:
                                            if kk > terminal_atom_arr[i]:
                                                chain_name = (
                                                    str(terminal_atom_arr[i])
                                                    + "-"
                                                    + str(kk)
                                                )
                                            else:
                                                chain_name = (
                                                    str(kk)
                                                    + "-"
                                                    + str(terminal_atom_arr[i])
                                                )
                                            if (
                                                chain_name
                                                not in length_between_terminal.keys()
                                            ):
                                                length_between_terminal[chain_name] = (
                                                    chain_atom_arr[:]
                                                )
                                                length_between_terminal[
                                                    chain_name
                                                ].append(kk)
                                            del_atom.append(kk)
                                            break_flag += 1
                                        else:
                                            chain_atom_arr.append(kk)
                                            j = kk
                                            break_flag_inner = 1
                                            break
                            if break_flag_inner == 1:
                                break
                            del_atom.append(chain_atom_arr[-1])
                            del chain_atom_arr[-1]
                    else:
                        j = k
                        chain_atom_arr.append(k)
                        break
            if break_flag == terminal_n - 1:
                break
    return length_between_terminal
